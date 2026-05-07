// tests/test_shared_het_rate.js

import { HET_RAMP, hetRateColor } from '../shared/het_rate.js';

let pass = 0, fail = 0;
function check(name, cond, detail = '') {
  if (cond) { console.log(`  ✓ ${name}${detail ? '  (' + detail + ')' : ''}`); pass++; }
  else      { console.log(`  ✗ ${name}  ${detail}`); fail++; }
}

console.log('--- HET_RAMP shape ---');
check('cold is rgb',                 HET_RAMP.cold.length === 3);
check('neutral is rgb',              HET_RAMP.neutral.length === 3);
check('warm is rgb',                 HET_RAMP.warm.length === 3);
check('cold = [33, 102, 172]',
      HET_RAMP.cold[0] === 33 && HET_RAMP.cold[1] === 102 && HET_RAMP.cold[2] === 172);
check('neutral = [247, 247, 247]',
      HET_RAMP.neutral[0] === 247 && HET_RAMP.neutral[1] === 247 && HET_RAMP.neutral[2] === 247);
check('warm = [178, 24, 43]',
      HET_RAMP.warm[0] === 178 && HET_RAMP.warm[1] === 24 && HET_RAMP.warm[2] === 43);
check('HET_RAMP frozen',             Object.isFrozen(HET_RAMP));

console.log('\n--- hetRateColor (sentinel cases) ---');
check('null → var(--ink-dimmer)',          hetRateColor(null)         === 'var(--ink-dimmer)');
check('undefined → var(--ink-dimmer)',     hetRateColor(undefined)    === 'var(--ink-dimmer)');
check('NaN → var(--ink-dimmer)',           hetRateColor(NaN)          === 'var(--ink-dimmer)');
check('Infinity → var(--ink-dimmer)',      hetRateColor(Infinity)     === 'var(--ink-dimmer)');

console.log('\n--- hetRateColor (anchor points) ---');
check('rate = 0   → cold',
      hetRateColor(0)   === `rgb(${HET_RAMP.cold[0]},${HET_RAMP.cold[1]},${HET_RAMP.cold[2]})`);
check('rate = 0.5 → neutral',
      hetRateColor(0.5) === `rgb(${HET_RAMP.neutral[0]},${HET_RAMP.neutral[1]},${HET_RAMP.neutral[2]})`);
check('rate = 1.0 → warm',
      hetRateColor(1.0) === `rgb(${HET_RAMP.warm[0]},${HET_RAMP.warm[1]},${HET_RAMP.warm[2]})`);

console.log('\n--- hetRateColor (clamps out of range) ---');
check('rate = -0.5 clamped to 0 (cold)',
      hetRateColor(-0.5) === `rgb(${HET_RAMP.cold[0]},${HET_RAMP.cold[1]},${HET_RAMP.cold[2]})`);
check('rate = 1.5 clamped to 1 (warm)',
      hetRateColor(1.5)  === `rgb(${HET_RAMP.warm[0]},${HET_RAMP.warm[1]},${HET_RAMP.warm[2]})`);

console.log('\n--- hetRateColor (interpolation midpoints) ---');
{
  // rate = 0.25: halfway cold→neutral. Each channel = (cold + neutral) / 2.
  const exp_r = Math.round((HET_RAMP.cold[0] + HET_RAMP.neutral[0]) / 2);
  const exp_g = Math.round((HET_RAMP.cold[1] + HET_RAMP.neutral[1]) / 2);
  const exp_b = Math.round((HET_RAMP.cold[2] + HET_RAMP.neutral[2]) / 2);
  check('rate = 0.25 → cold/neutral midpoint',
        hetRateColor(0.25) === `rgb(${exp_r},${exp_g},${exp_b})`);
}
{
  // rate = 0.75: halfway neutral→warm
  const exp_r = Math.round((HET_RAMP.neutral[0] + HET_RAMP.warm[0]) / 2);
  const exp_g = Math.round((HET_RAMP.neutral[1] + HET_RAMP.warm[1]) / 2);
  const exp_b = Math.round((HET_RAMP.neutral[2] + HET_RAMP.warm[2]) / 2);
  check('rate = 0.75 → neutral/warm midpoint',
        hetRateColor(0.75) === `rgb(${exp_r},${exp_g},${exp_b})`);
}

console.log('\n--- hetRateColor (output format) ---');
{
  // Sanity check on the format: matches rgb(N,N,N) with three integers in [0,255]
  const sample = hetRateColor(0.42);
  const m = sample.match(/^rgb\((\d+),(\d+),(\d+)\)$/);
  check('format = rgb(int,int,int)',  !!m, sample);
  if (m) {
    const [, r, g, b] = m.map(Number);
    check('all components 0-255',     [r,g,b].every(c => c >= 0 && c <= 255));
  }
}

console.log('\n=================');
console.log(`pass: ${pass}   fail: ${fail}`);
console.log('=================');
process.exit(fail === 0 ? 0 : 1);
