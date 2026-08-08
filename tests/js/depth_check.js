/* Depth must not be written as an inline opacity: an inline style beats the .node.dim class
   rule, which silently disabled every filter on the map. The entity count still changed, so
   it read as the toggles being broken rather than as a styling fault. */
const fs = require('fs');
const map = fs.readFileSync(process.argv[2], 'utf8');
const css = fs.readFileSync(process.argv[3], 'utf8');
const fail = [];
if (/style\.opacity\s*=/.test(map)) fail.push('map.js sets style.opacity directly');
if (!/setProperty\("--depth"/.test(map)) fail.push('map.js does not set --depth');
if (!/\.node\s*\{[^}]*opacity:\s*var\(--depth/.test(css)) fail.push('.node does not read --depth');
if (!/\.node\.dim\s*\{[^}]*opacity/.test(css)) fail.push('.node.dim no longer sets opacity');
if (fail.length) { console.error(fail.join('\n')); process.exit(1); }
console.log('ok');
