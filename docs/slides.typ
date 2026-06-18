#import "@preview/touying:0.7.4": *
#import themes.metropolis: *

#let ok(s) = text(fill: rgb("#2a7"), s)
#let bad(s) = text(fill: rgb("#c33"), s)
#let na(s) = text(fill: luma(45%), s)

// compact 3-col table: Circuit | qsyn | Vandaele(pyzx)
#let t-pyzx(..rows) = table(
  columns: (auto, auto, auto, auto),
  stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { (bottom: 0.2pt + luma(80%)) },
  inset: 4pt, align: (left, right, right, center),
  table.header(text(weight: "bold")[Circuit], text(weight: "bold")[old], text(weight: "bold")[new], text(weight: "bold")[Vandaele]),
  ..rows
)

// compact 4-col table: Circuit | paper | Vandaele | qsyn  (raw head-to-head)
#let t-raw(..rows) = table(
  columns: (auto, auto, auto, auto, auto),
  stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { (bottom: 0.2pt + luma(80%)) },
  inset: 4pt, align: (left, right, right, right, right),
  table.header(
    text(weight: "bold")[Circuit], text(weight: "bold")[paper],
    text(weight: "bold")[Vandaele], text(weight: "bold")[old], text(weight: "bold")[new]),
  ..rows
)

#let two(a, b) = grid(columns: (1fr, 1fr), column-gutter: 1.2em, a, b)

#show: metropolis-theme.with(
  aspect-ratio: "16-9",
  config-common(new-section-slide-fn: none),
  config-info(
    title: [Enhance FastTODD Algorithm in Qsyn],
    subtitle: [A 2 #sym.times 2 benchmark #sym.dash.em qsyn vs. the reference Rust implementation],
    author: [r14944052],
    date: datetime.today(),
    institution: [Quantum Circuit Optimization],
  ),
)

#title-slide()

== Setup

- *qsyn*: `tmerge #sym.arrow.r hopt #sym.arrow.r phasepoly fasttodd` (auto-decomposes Toffoli)
- *Vandaele*: `FastTMerge #sym.arrow.r InternalHOpt #sym.arrow.r FastTODD` (AVX2 build, paper-faithful)
- *Inputs*: pyzx-optimized (qsyn's set) and raw (Vandaele's set); 1800s + 12 GB caps
- Both single-threaded. Legend: #ok[T/s] done #sym.dot.c #bad[TLE/MLE/ERR] failed #sym.dot.c #na[N/A] unsupported input

= qsyn original (pyzx inputs)

== qsyn original #sym.dash.em Vandaele can't read these (1/2)

Vandaele rejects pyzx circuits (CZ gates) #sym.arrow.r #na[N/A] across the board.

#text(size: 9pt)[#two(
  t-pyzx(
    `adder_8`, ok[173/0s], ok[173/0s], na[N/A],
    `barenco_tof_3`, ok[16/0s], ok[16/0s], na[N/A],
    `barenco_tof_4`, ok[28/0s], ok[28/0s], na[N/A],
    `barenco_tof_5`, ok[40/0s], ok[40/0s], na[N/A],
    `barenco_tof_10`, ok[100/0s], ok[100/0s], na[N/A],
    `csla_mux_3`, ok[51/0s], ok[51/0s], na[N/A],
    `csum_mux_9`, ok[75/0s], ok[75/0s], na[N/A],
    `mod5_4`, ok[7/1s], ok[7/0s], na[N/A],
  ),
  t-pyzx(
    `mod_mult_55`, ok[35/0s], ok[35/0s], na[N/A],
    `mod_red_21`, ok[73/0s], ok[73/0s], na[N/A],
    `qcla_adder_10`, ok[155/0s], ok[156/0s], na[N/A],
    `qcla_com_7`, ok[95/0s], ok[95/0s], na[N/A],
    `qcla_mod_7`, ok[233/0s], ok[233/0s], na[N/A],
    `rc_adder_6`, ok[47/0s], ok[47/0s], na[N/A],
    `tof_3`, ok[15/0s], ok[15/0s], na[N/A],
  ),
)]

== qsyn original #sym.dash.em Vandaele can't read these (2/2)

#text(size: 9pt)[#two(
  t-pyzx(
    `tof_4`, ok[23/0s], ok[23/0s], na[N/A],
    `tof_5`, ok[31/0s], ok[31/0s], na[N/A],
    `tof_10`, ok[71/0s], ok[71/0s], na[N/A],
    `vbe_adder_3`, ok[24/0s], ok[24/0s], na[N/A],
    `gf2^4_mult`, ok[49/1s], ok[49/0s], na[N/A],
    `gf2^5_mult`, ok[77/1s], ok[77/0s], na[N/A],
    `gf2^6_mult`, ok[113/5s], ok[113/0s], na[N/A],
    `gf2^7_mult`, ok[149/17s], ok[143/2s], na[N/A],
  ),
  t-pyzx(
    `gf2^8_mult`, ok[203/69s], ok[203/3s], na[N/A],
    `gf2^9_mult`, ok[259/194s], ok[261/9s], na[N/A],
    `gf2^10_mult`, ok[293/486s], ok[283/41s], na[N/A],
    `gf2^16_mult`, bad[TLE], ok[791/1134s], na[N/A],
    `mod_adder_1024`, ok[1010/0s], ok[1010/0s], na[N/A],
    `nth_prime6`, ok[279/0s], ok[279/0s], na[N/A],
    `DEFAULT`, bad[TLE], bad[TLE], bad[MLE],
  ),
)]

= Vandaele paper (raw inputs)

== Vandaele paper #sym.dash.em head-to-head (1/2)

qsyn matches the paper's T; *our Vandaele beats it* (post-paper fix).

#text(size: 9pt)[#two(
  t-raw(
    `adder_8`, [170], ok[119/0s], ok[170/0s], ok[172/0s],
    `barenco_tof_3`, [16], ok[13/0s], ok[16/0s], ok[16/0s],
    `barenco_tof_4`, [28], ok[23/0s], ok[28/0s], ok[28/0s],
    `barenco_tof_5`, [40], ok[33/0s], ok[40/0s], ok[40/0s],
    `barenco_tof_10`, [100], ok[83/0s], ok[100/0s], ok[100/0s],
    `csla_mux_3`, [49], ok[39/0s], ok[53/0s], ok[53/0s],
    `csum_mux_9`, [73], ok[71/0s], ok[75/0s], ok[75/0s],
    `mod5_4`, [7], ok[7/0s], ok[7/0s], ok[7/0s],
  ),
  t-raw(
    `mod_mult_55`, [28], ok[17/0s], ok[32/0s], ok[32/0s],
    `mod_red_21`, [73], ok[51/0s], ok[73/0s], ok[73/0s],
    `qcla_adder_10`, [161], ok[109/2s], ok[156/0s], ok[160/0s],
    `qcla_com_7`, [95], ok[59/0s], ok[95/0s], ok[95/0s],
    `qcla_mod_7`, [237], ok[159/14s], ok[236/0s], ok[236/0s],
    `rc_adder_6`, [47], ok[37/0s], ok[47/0s], ok[47/0s],
    `tof_3`, [15], ok[13/0s], ok[15/0s], ok[15/0s],
  ),
)]

== Vandaele paper #sym.dash.em head-to-head (2/2)

`gf2` family: no paper number; qsyn is *~100#sym.times slower* and TLEs `gf2^16`.

#text(size: 9pt)[#two(
  t-raw(
    `tof_4`, [23], ok[19/0s], ok[23/0s], ok[23/0s],
    `tof_5`, [31], ok[25/0s], ok[31/0s], ok[31/0s],
    `tof_10`, [71], ok[55/0s], ok[71/0s], ok[71/0s],
    `vbe_adder_3`, [24], ok[19/0s], ok[24/0s], ok[24/0s],
    `gf2^4_mult`, na[#sym.dash.en], ok[49/0s], ok[49/1s], ok[49/0s],
    `gf2^5_mult`, na[#sym.dash.en], ok[81/0s], ok[75/1s], ok[75/0s],
    `gf2^6_mult`, na[#sym.dash.en], ok[113/0s], ok[115/5s], ok[105/0s],
  ),
  t-raw(
    `gf2^7_mult`, na[#sym.dash.en], ok[155/1s], ok[147/21s], ok[147/1s],
    `gf2^8_mult`, na[#sym.dash.en], ok[205/1s], ok[201/41s], ok[201/4s],
    `gf2^9_mult`, na[#sym.dash.en], ok[257/2s], ok[255/197s], ok[243/11s],
    `gf2^10_mult`, na[#sym.dash.en], ok[315/5s], ok[315/205s], ok[319/13s],
    `gf2^16_mult`, na[#sym.dash.en], ok[797/172s], bad[TLE], ok[807/1075s],
    `mod_adder_1024`, [1009], bad[TLE], ok[1009/1s], ok[1009/1s],
    `DEFAULT`, [39666], bad[MLE], bad[TLE], bad[TLE],
  ),
)]

= Vandaele new (added post-paper)

== Vandaele new #sym.dash.em the gadgetization wall

Hadamard-dense circuits: gadgetize-all makes Vandaele *TLE/MLE*; qsyn finishes (shallow).

#text(size: 9pt)[#two(
  t-raw(
    `grover_5`, [166], ok[143/2s], ok[166/0s], ok[166/0s],
    `ham15-low`, [94], ok[77/0s], ok[97/0s], ok[97/0s],
    `ham15-med`, [212], ok[137/8s], ok[212/0s], ok[212/0s],
    `hwb6`, [75], ok[51/0s], ok[72/0s], ok[72/0s],
    `qft_4`, [66], ok[53/0s], bad[ERR], bad[ERR],
  ),
  t-raw(
    `ham15-high`, [1019], bad[TLE], ok[1011/0s], ok[1011/0s],
    `hwb8`, [3487], bad[TLE], ok[3502/1s], ok[3505/1s],
    `hwb10`, [15693], bad[MLE], ok[15816/11s], ok[15840/10s],
    `hwb11`, [44188], bad[MLE], ok[44356/54s], ok[44417/51s],
  ),
)]

= Findings

== Findings

- *Input compatibility*: Vandaele can't ingest qsyn's pyzx circuits (CZ gates) #sym.dash.em #na[N/A] on all 47.
- *Speed (the fix)*: old qsyn was *~70#sym.dash.en#h(0pt)100#sym.times slower* on `gf2` and TLE'd `gf2^16`. Ported FastTODD (`new` col) cuts this *~15#sym.dash.en#h(0pt)50#sym.times*: gf2^9 197#sym.arrow.r#h(0pt)11s, gf2^10 205#sym.arrow.r#h(0pt)13s. `gf2^16` is *borderline* (this sweep 1075#sym.dash.en#h(0pt)1134s; an earlier run hit 1861s = TLE) #sym.dash.em sits right at the 1800s cap, not robustly fixed, still ~6#sym.times Vandaele's 172s.
- *T-count #sym.tilde unchanged*: ported tracks old to within a handful (adder_8 170#sym.arrow.r#h(0pt)172, gf2^10 pyzx 293#sym.arrow.r#h(0pt)283); the speedup costs ~no quality.
- *Gadgetization*: Hadamard-dense circuits (hwb, ham15) make Vandaele *TLE/MLE* #sym.dash.em but these are post-paper inputs, not the paper's setup.
- *Quality*: our HEAD Vandaele *beats* the paper's FastTODD T-counts; qsyn matches the paper.
- *`DEFAULT` (640q)*: defeats both #sym.dash.em qsyn TLE, Vandaele MLE (>52 GB).
- *Failure modes are opposite*: qsyn = time-bound (TLE, ~5.5 GB); Vandaele = memory-bound (MLE).
