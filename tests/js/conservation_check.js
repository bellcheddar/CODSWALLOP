/* The entity-index-to-seed-position frame used to colour a structure by conservation.

   This is the same conversion that shipped wrong in the contacts artefact for 52 of 71
   families and looked entirely plausible on screen, so it is scored on residue identity
   rather than assumed from query_beg, and refused when it does not agree. */
const fs = require('fs');
const src = fs.readFileSync(process.argv[2], 'utf8');
const fail = [];

function lift(name) {
  const i = src.indexOf('function ' + name + '(');
  if (i < 0) { fail.push('viewer.js no longer defines ' + name); return 'function(){}'; }
  let depth = 0, j = src.indexOf('{', i);
  for (; j < src.length; j++) {
    if (src[j] === '{') depth++;
    else if (src[j] === '}' && --depth === 0) break;
  }
  return src.slice(i, j + 1);
}
function num(name) {
  const m = new RegExp('var\\s+' + name + '\\s*=\\s*([0-9.]+)').exec(src);
  if (!m) { fail.push('viewer.js no longer defines ' + name); return null; }
  return Number(m[1]);
}

global.OFFSET_WINDOW = num('OFFSET_WINDOW');
const MIN_AGREEMENT = num('MIN_AGREEMENT');
const bestOffset = eval('(' + lift('bestOffset') + ')');

const SEED = 'MSHHWGYGKHNGPEHWHKDFPIAKGERQSPVDIDTHTAKYDPSLKPLSVSYDQATSLRILNNGHAFNVEFDDSQDKAVLK'
           + 'GGPLDGTYRLIQFHFHWGSLDGQGSEHTVDKKKYAAELHLVHWNTKYGDFGKAVQQPDGLAVLGIFLKVGSAKPGLQKVV';
const columns = SEED.split('').map((aa, i) => ({ pos: i + 1, seed: aa, conservation: 0.5 }));

function chain(seq, labelStart) {
  return seq.split('').map((aa, i) => ({ chain: 'A', auth: labelStart + i,
                                         label: labelStart + i, aa }));
}
function check(c, m) { if (!c) fail.push(m); }

// --- the construct starts where the alignment says: offset 0 -----------------------------
{
  const fit = bestOffset(chain(SEED, 1), columns, 0);
  check(fit.offset === 0, 'expected offset 0, got ' + fit.offset);
  check(fit.score === 1, 'a perfect match should score 1, got ' + fit.score);
}

// --- an N-terminal tag shifts the entity numbering away from the expected offset ----------
// This is why the offset is searched rather than taken from query_beg: the expected value is
// wrong by exactly the length of the tag, and a linear map would frame-shift the whole
// protein while still producing a plausible-looking coloured picture.
{
  const TAG = 'MGSSHHHHHHSSGLVPRGSH';
  const fit = bestOffset(chain(TAG + SEED, 1), columns, 0);
  check(fit.offset === -TAG.length,
        'expected the tag to be absorbed as offset ' + (-TAG.length) + ', got ' + fit.offset);
  check(fit.score > 0.97, 'the tagged construct should still match well, got ' + fit.score);
}

// --- a construct numbered from 1 on a fragment of the seed --------------------------------
{
  const frag = SEED.slice(40, 120);
  const fit = bestOffset(chain(frag, 1), columns, 0);
  check(fit.offset === 40, 'expected offset 40, got ' + fit.offset);
  check(fit.score > 0.97, 'the fragment should match well, got ' + fit.score);
}

// --- an unrelated protein must NOT be given a frame ---------------------------------------
// The whole point of scoring on identity: without the agreement floor this returns whichever
// offset happened to score least badly, and the structure gets coloured with somebody else's
// conservation.
{
  let other = '';
  for (let i = 0; i < 120; i++) other += 'AAAGGGSSS'[i % 9];
  const fit = bestOffset(chain(other, 1), columns, 0);
  check(fit.score < MIN_AGREEMENT,
        'an unrelated chain scored ' + fit.score + ', at or above the ' + MIN_AGREEMENT +
        ' floor, so it would have been coloured');
}

// --- a handful of point mutations must not lose the frame ---------------------------------
{
  const mutated = SEED.split('');
  [10, 50, 90, 130].forEach(i => { mutated[i] = mutated[i] === 'A' ? 'W' : 'A'; });
  const fit = bestOffset(chain(mutated.join(''), 1), columns, 0);
  check(fit.offset === 0, 'four mutations should not move the frame, got ' + fit.offset);
  check(fit.score > MIN_AGREEMENT, 'a mutant should still be coloured, got ' + fit.score);
}

// --- the renderer must actually refuse a bad frame -----------------------------------------
check(/fit\.score\s*<\s*MIN_AGREEMENT/.test(src),
      'conservationColour no longer refuses a low-agreement frame');
check(/col\.seed\s*!==\s*residues\[r\]\.aa/.test(src),
      'residues that differ from the seed are no longer left uncoloured');

if (fail.length) { console.error(fail.join('\n')); process.exit(1); }
console.log('ok');
