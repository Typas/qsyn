#import "@preview/fletcher:0.5.8"
#import "@preview/tblr:0.5.0"

#set text(size: 12pt,
    font: (
    "STIX Two Text",
    "思源宋體",
))
#show math.equation: set text(size: 12pt,
    font: (
    "STIX Two Math",
    "思源宋體",
))

#set document(
  title: [Final Project: FastTODD Implementation],
  author: "R14944052 Chen-Hao Liao",
  date: datetime.today(),
)

#align(center)[
  #title()
  R14944052 Chen-Hao Liao #h(1.5em) #datetime.today().display("[year]-[month]-[day]")
]

// Render a benchmark CSV as a booktabs-style tblr table.
#let bench(path, cols, aligns, cap) = tblr.tblr(
  columns: cols,
  header-rows: 1,
  align: aligns,
  caption: cap,
  placement: none,
  tblr.rows(within: "header", 0, hooks: strong),
  content-hook: tblr.from-csv,
  read(path),
)

= Introduction

#bibliography("references.bib", style: "ieee")
