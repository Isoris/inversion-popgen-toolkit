// build_synth_kmeans.mjs — synthesize a 226-sample, 3-window, K=3 kmeans
// labels JSON matching the LG28 prototype's karyotype shape (60/106/60).
// Cluster ids are window-local (as in real K-means output).
// Window 0: cluster 0 = HOM_REF (60), cluster 1 = HET (106), cluster 2 = HOM_INV (60)
// Window 1: cluster 1 = HOM_REF (60), cluster 2 = HET (106), cluster 0 = HOM_INV (60)
// Window 2: cluster 2 = HOM_REF (60), cluster 0 = HET (106), cluster 1 = HOM_INV (60)
// (intentionally permute cluster ids across windows — that's the realistic
//  case where K-means doesn't preserve label identity)

import fs from 'node:fs';

const N_SAMPLES = 226;
const HOM_REF = 60, HET = 106, HOM_INV = 60;
console.assert(HOM_REF + HET + HOM_INV === N_SAMPLES);

// True karyotype per sample (for our reference, not in the file):
//   samples 0..59      = HOM_REF
//   samples 60..165    = HET
//   samples 166..225   = HOM_INV
const truth = [];
for (let i = 0; i < N_SAMPLES; i++) {
  if (i < HOM_REF) truth.push('HOM_REF');
  else if (i < HOM_REF + HET) truth.push('HET');
  else truth.push('HOM_INV');
}

// Cluster-id permutation per window
const perms = [
  { HOM_REF: 0, HET: 1, HOM_INV: 2 },   // window 0
  { HOM_REF: 1, HET: 2, HOM_INV: 0 },   // window 1
  { HOM_REF: 2, HET: 0, HOM_INV: 1 },   // window 2
];

const windows = perms.map((perm, w) => {
  const labels = truth.map(t => perm[t]);
  return {
    w,
    start_bp: 15115000 + w * 1000000,
    end_bp:   15115000 + (w + 1) * 1000000,
    band_labels: labels,
  };
});

const out = {
  _comment: 'Synthetic LG28 prototype: 226 samples, 60/106/60 karyotype, 3 windows, K=3 per window with permuted cluster ids.',
  n_samples: N_SAMPLES,
  K_locus: 9,                    // 3 windows × 3 clusters/window
  windows,
};
fs.writeFileSync('synth_LG28_kmeans.json', JSON.stringify(out, null, 2));
console.log('wrote synth_LG28_kmeans.json');
