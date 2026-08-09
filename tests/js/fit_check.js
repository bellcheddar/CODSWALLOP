/* The panel fit, and the clamp that keeps an outlier from setting the scale.

   The case this exists for: hen lysozyme's cryo-EM amyloid fibrils (9J0L, 9J0M) are the same
   sequence as the seed and TM 0.15-0.20 to all 76 other representatives. Fitting the panel to
   them put 1,686 of 1,688 nodes inside 5.4 % of its width.

   map.js is a browser IIFE with no exports, so the two pure functions are lifted out by
   source and evaluated. That is deliberate: reimplementing them here would test the copy. */
const fs = require('fs');
const src = fs.readFileSync(process.argv[2], 'utf8');
const fail = [];

function lift(name) {
  const i = src.indexOf('function ' + name + '(');
  if (i < 0) { fail.push('map.js no longer defines ' + name); return null; }
  // Balance braces from the function's opening brace to its close.
  let depth = 0, start = src.indexOf('{', i), j = start;
  for (; j < src.length; j++) {
    if (src[j] === '{') depth++;
    else if (src[j] === '}' && --depth === 0) break;
  }
  return src.slice(i, j + 1);
}

function constant(name) {
  const m = new RegExp('var\\s+' + name + '\\s*=\\s*([0-9.]+)').exec(src);
  if (!m) { fail.push('map.js no longer defines ' + name); return null; }
  return Number(m[1]);
}

const robustExtent = eval('(' + lift('robustExtent') + ')');
const FIT_QUANTILE = constant('FIT_QUANTILE');
const FIT_HEADROOM = constant('FIT_HEADROOM');
const FIT_MIN_GAIN = constant('FIT_MIN_GAIN');
const FIT_OVERFLOW = constant('FIT_OVERFLOW');
// The lifted function closes over these, so they have to exist in this scope by those names.
global.FIT_QUANTILE = FIT_QUANTILE;
global.FIT_HEADROOM = FIT_HEADROOM;
global.FIT_MIN_GAIN = FIT_MIN_GAIN;

function check(cond, msg) { if (!cond) fail.push(msg); }

// --- a well-behaved family is fitted to its true extent, and pins nothing ---------------
{
  const vals = [];
  for (let i = 0; i < 80; i++) vals.push(0.2 + 0.8 * (i / 79));   // a smooth spread to 1.0
  const ext = robustExtent(vals);
  check(Math.abs(ext - 1.0) < 1e-9,
    'an even spread should fit to its maximum, got ' + ext);
  check(!vals.some(v => v / ext > FIT_OVERFLOW),
    'an even spread must not pin anything');
}

// --- the lysozyme shape: 75 points in a tight cloud, one 15x further out ----------------
{
  const vals = [];
  for (let i = 0; i < 75; i++) vals.push(0.02 + 0.045 * (i / 74));  // cloud out to 0.065
  vals.push(1.0);                                                    // the amyloid fibrils
  const ext = robustExtent(vals);
  check(ext < 0.2, 'the outlier still set the scale: extent ' + ext);
  check(1.0 / ext > 5, 'the outlier should be well beyond the fitted extent');
  const pinned = vals.filter(v => v / ext > FIT_OVERFLOW).length;
  check(pinned === 1, 'exactly the outlier should be pinned, got ' + pinned);
  // The bulk has to actually gain the panel back.
  check(0.065 / ext > 0.4,
    'the cloud should fill most of the fitted range, got ' + (0.065 / ext).toFixed(2));
}

// --- a long tail is a tail, not an outlier: nothing is pinned ---------------------------
{
  const vals = [];
  for (let i = 0; i < 60; i++) vals.push(0.1 + 0.5 * (i / 59));
  for (let i = 0; i < 20; i++) vals.push(0.6 + 0.4 * (i / 19));   // 25 % of it in the tail
  const ext = robustExtent(vals);
  const pinned = vals.filter(v => v / ext > FIT_OVERFLOW).length;
  check(pinned === 0, 'a long tail must not be treated as outliers, pinned ' + pinned);
}

// --- degenerate inputs must not divide by zero ------------------------------------------
check(robustExtent([]) === 0, 'an empty family should give a zero extent, not NaN');
check(robustExtent([0, 0, 0]) === 0, 'every node stacked on one point should give zero');
check(robustExtent([0.5]) === 0.5, 'a single point fits to itself');

// --- the clamp keeps its bearing and lands exactly on the rim ----------------------------
{
  const frame = { cx: 100, cy: 50, scale: 1, limX: 40, limY: 20 };
  const place = function (project, n) {
    const p = project(n);
    let dx = p.x - frame.cx, dy = p.y - frame.cy;
    const over = Math.max(Math.abs(dx) / frame.limX, Math.abs(dy) / frame.limY);
    if (over > 1) { dx /= over; dy /= over; }
    return { x: frame.cx + dx, y: frame.cy + dy, z: p.z };
  };
  const id = n => ({ x: n.x, y: n.y, z: 0 });

  const inside = place(id, { x: 120, y: 55 });
  check(inside.x === 120 && inside.y === 55, 'a point inside the panel must not be moved');

  const out = place(id, { x: 500, y: 60 });                 // far right, slightly down
  check(Math.abs(out.x - (frame.cx + frame.limX)) < 1e-9,
    'a clamped point should land on the rim, got x=' + out.x);
  const bearingIn = Math.atan2(60 - frame.cy, 500 - frame.cx);
  const bearingOut = Math.atan2(out.y - frame.cy, out.x - frame.cx);
  check(Math.abs(bearingIn - bearingOut) < 1e-9,
    'clamping must preserve the bearing: ' + bearingIn + ' vs ' + bearingOut);

  // Both axes over at once: whichever is worse decides, and neither ends up outside.
  const corner = place(id, { x: 100 + 400, y: 50 + 400 });
  check(Math.abs(corner.x - frame.cx) <= frame.limX + 1e-9
     && Math.abs(corner.y - frame.cy) <= frame.limY + 1e-9,
    'a diagonal overflow escaped the panel');
}

// --- the renderer must actually use it ---------------------------------------------------
check(/Constellation\.prototype\.place\s*=/.test(src), 'place() is gone from map.js');
check(!/var\s+p3\s*=\s*project\(n\)/.test(src),
  'the node loop projects without clamping');
check(/robustExtent\(toDraw/.test(src),
  'the fit no longer measures the drawn points');
check(/offscalemark/.test(src), 'the off-scale chevron is gone');

if (fail.length) { console.error(fail.join('\n')); process.exit(1); }
console.log('ok');
