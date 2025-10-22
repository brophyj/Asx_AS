# Install if needed

load_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

load_package("DiagrammeRsvg")
load_package("DiagrammeR")
load_package("rsvg")



# Create the graph and save as object
graph <- grViz("
digraph prisma_flow {
  graph [
    layout = dot,
    rankdir = TB,
    labelloc = \"t\",  // places the title at the top
    label = \"Figure 1 PRISMA Flowchart of Study Identification and Inclusion\",
    fontsize = 18,
    fontname = \"Helvetica-Bold\"
  ]
  node [shape = box, style = rounded, fontsize = 10, fontname = Helvetica]

  A [label = '4 RCTs identifed from\\n 4 previous systematic review \\n until 2024/4/15
   \\nDOI {10.1016/j.amjcard.2025.02.025} \\nDOI: 10.1016/j.jacc.2024.11.006  \\nDOI: {10.1016/j.amjcard.2025.02.025} \\nDOI: {10.1016/j.amjcard.2025.02.025}
']
  B [label = 'Additional records identified\\nfrom 2024/4/15 until 2025/10/15\\n(n=93)']
  C [label = 'Records screened\\n(n = 97)']
  D [label = 'Records excluded\\n(n = 89)']
  E [label = 'Full-text articles assessed\\nfor eligibility\\n(n = 4)']
  F [label = 'Studies included in\\nquantitative synthesis\\n(meta-analysis)\\n(n = 4)']
  G [label = 'Search strategy:\\nasymptomatic AND severe AND \naortic stenosis AND \n(aortic valve replacement \nOR SAVR OR TAVR OR TAVI OR \ntranscatheter aortic valve OR \nconservative management OR \nwatchful waiting ']

  A -> C
  B -> C
  C -> D
  C -> E
  E -> F
  G
}
")

# Convert to SVG and then to PNG and PDF
svg <- export_svg(graph)
charToRaw(svg) %>% rsvg_png("Output/prisma_flowchart.png", width = 1200)
charToRaw(svg) %>% rsvg_pdf("Output/prisma_flowchart.pdf")

