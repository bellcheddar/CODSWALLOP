/* Superposition geometry, checked against an answer known in advance.
 *
 * Both bugs this catches produced a picture that looked like a superposition:
 *   - testing the handedness of U instead of det(U.V^T) returned a reflection (det -1),
 *   - building U.V^T instead of V.U^T returned the rotation transposed.
 * Each gave 21 A RMSD on two clouds related by an exact rotation, and neither would have
 * been obvious on screen: the model lands near the structure and looks roughly aligned.
 */
const fs = require('fs');
const src = fs.readFileSync(process.argv[2], 'utf8');

function grab(name) {
  const i = src.indexOf('function ' + name + '(');
  if (i < 0) throw new Error('missing function ' + name);
  let d = 0, end = i;
  for (let k = src.indexOf('{', i); k < src.length; k++) {
    if (src[k] === '{') d++;
    else if (src[k] === '}') { d--; if (d === 0) { end = k; break; } }
  }
  return src.slice(i, end + 1);
}
eval(grab('kabsch') + grab('jacobi') + grab('dot') + grab('cross') + grab('det3'));

function det(m) {
  return m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
       - m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
       + m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
}
const fail = [];

// 1. Two clouds related by an exact rotation must fit with zero error.
for (const deg of [15, 40, 90, 175]) {
  const th = deg * Math.PI / 180, c = Math.cos(th), s = Math.sin(th);
  const R = [[c, -s, 0], [s, c, 0], [0, 0, 1]], T = [7, -3, 11];
  const mob = Array.from({length: 80}, (_, i) =>
    [Math.sin(i) * 20, Math.cos(i * 1.7) * 18, (i % 13) * 3]);
  const ref = mob.map(p => [0, 1, 2].map(k =>
    R[k][0]*p[0] + R[k][1]*p[1] + R[k][2]*p[2] + T[k]));
  const fit = kabsch(mob, ref);
  if (!fit || fit.rmsd > 1e-6) fail.push(`${deg} deg: rmsd ${fit && fit.rmsd}`);
  if (!fit || Math.abs(det(fit.u) - 1) > 1e-6) fail.push(`${deg} deg: det ${fit && det(fit.u)}`);
  if (fit) {
    for (let k = 0; k < 3; k++) {
      if (Math.abs(fit.t[k] - T[k]) > 1e-6) fail.push(`${deg} deg: t[${k}] ${fit.t[k]}`);
      for (let m = 0; m < 3; m++) {
        if (Math.abs(fit.u[k][m] - R[k][m]) > 1e-6) fail.push(`${deg} deg: u[${k}][${m}]`);
      }
    }
  }
}

// 2. A mirrored cloud must NOT be fitted by a reflection: the best rotation is a poor fit,
//    and returning det -1 would make it look perfect.
const base = Array.from({length: 60}, (_, i) => [Math.sin(i)*20, Math.cos(i*1.3)*15, i]);
const mirrored = base.map(p => [p[0], p[1], -p[2]]);
const m = kabsch(base, mirrored);
if (!m || Math.abs(det(m.u) - 1) > 1e-6) fail.push('mirror: det is not +1');
if (m && m.rmsd < 1) fail.push('mirror: fitted a reflection as though it were a rotation');

// 3. Fewer than three points has no unique answer.
if (kabsch([[0,0,0],[1,0,0]], [[0,0,0],[1,0,0]]) !== null) fail.push('two points returned a fit');

if (fail.length) { console.error(fail.join('\n')); process.exit(1); }
console.log('ok');
