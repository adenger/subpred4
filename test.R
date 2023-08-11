library(GOSemSim)
library('org.Hs.eg.db')
d <- godata('org.Hs.eg.db', ont="MF", computeIC=FALSE)  
goSim("GO:0000006", "GO:0000007", semData=d, measure="Wang")