
library(ggplot2)
library(ggtree)
library(phytools)
library(ape)


tree_bl <- read.tree("../tree_pnas_mech_nucl.txt")
plot(tree_bl)
tiplabels()

tree_bl$tip.label

keep_species <- c("Mon_2012_244_Mech-messenoides-messenoides_Ec-WA", 
 "CAM071976_Mech-polymnia-chimborazona_Ec-WEST",
 "CAM028331_Mech-polymnia-dorissides_Ec-WA",
 "Sal_1416_Mech-polymnia-dorissides_Col-WA",
 "CAM028192_Mech-lysimnia-roqueensis_Ec-WA",
 "CAM028198_Mech-lysimnia-roqueensis_Ec-WA",
 "CAM008476_Mech-macrinus-macrinus_Pan-WEST",
 "Mon_2012_197_Mech-mazaeus-fallax_Ec-WA",
 "CAM070759_Mech-menapis-mantineus_Ec-WEST",
"CAM070479_For-proceris_Ec-WA")

tree_pruned <- keep.tip(tree_bl, keep_species)
plot(tree_pruned)
nodelabels()
tiplabels()

tree_pruned <- reroot(tree_pruned, node.number = 7 , position = 0.02)
plot(tree_pruned)
nodelabels()

ggtree(tree_pruned) +
  geom_nodelab() +
  geom_tiplab() +
  xlim(c(0,0.1))


tree_pruned$tip.label <- c("ilMecMess1", "ilMecPoly75","ilMecPoly1", "ilMecPoly126", "ilMecLysi317", "ilMecLysi212", "ilDanPlex1",  "ilMecMacr1",   "ilMecMaza1","ilMecMena1" )

write.tree(tree_pruned, file = "tree_pnas_pruned.txt")

tree_pruned$tip.label
[1] "Mon_2012_244_Mech-messenoides-messenoides_Ec-WA"
[2] "CAM071976_Mech-polymnia-chimborazona_Ec-WEST"   
[3] "CAM028331_Mech-polymnia-dorissides_Ec-WA"       
[4] "Sal_1416_Mech-polymnia-dorissides_Col-WA"       
[5] "CAM028192_Mech-lysimnia-roqueensis_Ec-WA"       
[6] "CAM070479_For-proceris_Ec-WA"                   
[7] "CAM008476_Mech-macrinus-macrinus_Pan-WEST"      
[8] "Mon_2012_197_Mech-mazaeus-fallax_Ec-WA"         
[9] "CAM070759_Mech-menapis-mantineus_Ec-WEST"


