
# ---- AI ----
# print(col_ext_full[reg=="ai", list(spp, orig_tax_id, full_firstYear, full_lastYear, full_propYear)], nrow=Inf)
man_remove_ai <- c(
"Aphrocallistes vastus", # missing until 1994, 5%, is sponge; # ---- TOSS ----
"Aplidium soldatovi", # missing until 2006, shows up at 15%, is a sea squirt # ---- TOSS ----
# "Aptocyclus ventricosus", # missing only first year, 10%, fish (lump sucker), --- KEEP
# "Bathyraja panthera", # 2010, high certainty, is skate # ---- KEEP ----
# "Bathyraja taranetzi", # 1991, 5%, high certainty, skate, # ---- KEEP ---- ???
# "Careproctus rastrinus", # present first year, high certainty all years, is snailfish, # ---- KEEP ----
"Ceramaster japonicus", # first in 1994 (1991 full), 5%, is sea star, uncertain all yrs; # ---- TOSS ----
"Ceramaster patagonicus", # first 1991, 15%, sea tar, uncertain all yrs; # ---- TOSS ----
# "Cheiraster dawsoni", # first 1991, 5%, sea star, certain all yrs; # ---- KEEP ----
"Chlamys albida", # first 1997, 15%, scallop, uncertain all yrs (4); # ---- TOSS ----
# "Chlamys rubida", # first 1986, 5%, scallop, uncertain all yrs, # ---- KEEP ----
"Chrysaora melanaster", # shows up 2000, 5%, nettle, high certainty # ---- TOSS ---- ????
"Clathrina blanca", # shows up 1997, 30%, sponge, low certainty all yrs; # ---- TOSS ----
"Crossaster papposus", # shows up 1991 ('87 goa), 10%, high certainy all yrs; # ---- TOSS ----
"Cucumaria fallax", # 1991, 10%, cucumber, high certainty all yrs; probably orig ID as Holothuroidea; # ---- TOSS ----
"Diplopteraster multipes", # 1994 ('91 full), 25%, sea star, high certainty; # ---- TOSS ----
"Elassochirus cavimanus", # 1997, 15%, hermit, high certainty; # ---- TOSS ----
# "Enteroctopus dofleini", # 1986 (1980 [1st yr] full), giant octo, high certainty all yrs; # ---- KEEP ----
"Fanellia compressa", # 1996, 20%, sea fan, low certainty all yrs; # ---- TOSS ----
# "Fusitriton oregonensis", # 1986 (1980 full), snail, 40%, high certainty all yr; # ---- KEEP ----
"Geodia lendenfeldi", # 2004, 10%, sponge, low certainty; # ---- TOSS ----
"Halichondria panicea", # 1994, 25%, sponge, low certainty; # ---- TOSS ----
"Halichondria sitiens", # 1997, 50%, sponge, low certainty; # ---- TOSS ----
"Halocynthia aurantium", # 1991, 5%, sea peach, high certainty; # ---- TOSS ----
# "Hemilepidotus zapus", # 1983, 5%, fish (irish lord), high certainty; # ---- KEEP ----
"Henricia aspera", # 1997, 3%, sea star, low certainty; # ---- TOSS ----
"Henricia leviuscula", # 1991, sea star, 1%, low certainty; # ---- TOSS ----
"Hippasteria phrygiana", # 1991, sea star, 10%, recorded as 2 spp; low certainty; # ---- TOSS ----
"Homaxinella amphispicula", # 1997, 25%, firm finger sponge, low certainty; # ---- TOSS ----
# "Hyas lyratus", # 1986 (1980 full), 5%, crab, high certainty; # ---- KEEP ----
"Isodictya rigida", # 1997, 10%, sponge, low certainty; # ---- TOSS ----
"Lebbeus groenlandicus", # 1994, 5%, shrimp, high certainty --- shrimp UNID until 1993; # ---- TOSS ----
"Leptasterias coei", # 2002, 5%, sea star, 2 spp, low certainty # ---- TOSS ----
# "Lethotremus muticus", # 1983 (1980 full), snailfish, medium certainty; # ---- KEEP ----
# "Modiolus modiolus", # 1986, mussel, 5%, high certainty; # ---- KEEP ---- ?????
"Muriceides nigra", # 2004, soft coral, 5%, low certainty; # ---- TOSS ----
"Mycale loveni", # 1991, tree sponge, 5%, low certainty; # ---- TOSS ----
"Ophiura sarsii", # 1994, sea star, 5%, high certainty; # ---- TOSS ----
# "Pagurus aleuticus", # 1991 aleutian hermit; 5%, high certainty; # ---- KEEP ---- ???
"Pagurus brandti", # 1991 sponge hermit, 5%, low certainty; # ---- TOSS ----
"Pagurus confragosus", # 2000, hermit, 5%, medium certainty; # ---- TOSS ----
"Pagurus trigonocheirus", # 2000, hermit, 5%, low certainty; # ---- TOSS ----
# "Pandalus tridens", # 1994 (1980 full), shrimp, 10%, medium certainty; # ---- KEEP ---- ???
"Plakina tanaga", # 2006, sponge, 5%, medium certainty; # ---- TOSS ----
"Pododesmus macrochisma", # 1986, alaska jingle, 5%, 2 spp, high certainty; # ---- TOSS ---- ???
"Porella compressa", # 1994, sea fan, 5%, low certainty; # ---- TOSS ----
"Pseudarchaster parelii", # 1994, sea star, 15%, low certainty; # ---- TOSS ----
"Pteraster marsippus", # 1997, sea star, 15%, low certainty; # ---- TOSS ----
"Pteraster tesselatus", # 1991, 10%, sea star, low certainty; # ---- TOSS ----
"Sebastes melanostictus", # 2006, 25%, rockfish, previous grouped w/ S. aleutians # ---- TOSS ----
"Stegophiura ponderosa", # 1994, 10%, sea star, high certainty; # ---- TOSS ----
# "Strongylocentrotus droebachiensis", # leaves 2012, present all yrs in full, urchin, low certainty # ---- KEEP ----
"Strongylocentrotus polyacanthus", # 2002 40% urchin, low certainty; # ---- TOSS ----
"Styela rustica", # 1994, sea potato, 15%, high certainty; # ---- TOSS ----
"Suberites ficus", # 1991, sea sponge, 15%, high certainty; # ---- TOSS ----
# "Triglops macellus", # present 1st yr, abs 1986, very rare (~3%), high certainty; # ---- KEEP ----
""
)


# ---- EBS ----
# print(col_ext_full[reg=="ebs", list(spp, orig_tax_id, full_firstYear, full_lastYear, full_propYear)], nrow=Inf)
man_remove_ebs <- c(
# "Aphrodita negligens",
# "Argis dentata",
# "Argis lar",
# "Artediellus pacificus",
# "Aurelia labiata",
# "Buccinum angulosum",
# "Buccinum polare",
# "Careproctus phasma",
# "Chrysaora melanaster",
# "Ciliatocardium ciliatum",
# "Colus herendeenii",
# "Crangon dalli",
# "Cryptonatica russa",
# "Cyanea capillata",
# "Echinarachnius parma",
# "Eleginus gracilis",
# "Enteroctopus dofleini",
# "Eunoe depressa",
# "Euspira pallida",
# "Gersemia rubiformis",
# "Glebocarcinus oregonensis",
# "Grandicrepidula grandis",
# "Halichondria panicea",
# "Halocynthia aurantium",
# "Hiatella arctica",
# "Icelus spatula",
# "Leptasterias groenlandica",
# "Liparis gibbus",
# "Liponema brevicorne",
# "Mactromeris polynyma",
# "Metridium dianthus",
# "Metridium farcimen",
# "Musculus discors",
# "Neocrangon communis",
# "Pyrulofusus melonis",
# "Rhamphostomella costata",
# "Serratiflustra serrulata",
# "Serripes notabilis",
# "Siliqua alta",
# "Stomphia coccinea",
# "Styela rustica",
# "Tellina lutea",
# "Trichodon trichodon",
# "Triglops pingelii",
# "Triglops scepticus",
# "Tritonia diomedea",
# "Urticina crassicornis",
# "Volutopsius fragilis",
# "Volutopsius middendorffii",
""
)

# ---- GMex ----
# print(col_ext_full[reg=="gmex", list(spp, orig_tax_id, full_firstYear, full_lastYear, full_propYear)], nrow=Inf)
man_remove_gmex <- c(
# "Anadara baughmani", # 1985, comes and goes, 5%, then booms to 10%, ark, # ---- KEEP ----
# "Anchoviella perfasciata", # 1985, anchovy, 10%, # ---- KEEP ----
# "Arenaeus cribrarius", # 1985, crab, 5%, # ---- KEEP ----
"Astropecten cingulatus", # 1993, sea star, 15%, # ---- TOSS ----
"Astropecten duplicatus", # 1992, sea star, 15% # ---- TOSS ----
"Aurelia aurita", # 1993, moon jelly, rester says ID'd, but rel abund bad; # ---- TOSS ----
# "Chaetodipterus faber", # pres 1st yr, spadefish, # ---- KEEP ----
# "Chrysaora quinquecirrha", # 1989, 5%, nettle, # ---- KEEP ----
# "Citharichthys macrops", # pres 1st yr, flatfish, # ---- KEEP ----
# "Decapterus punctatus", # pres 1st yr, scad, # ---- KEEP ----
# "Engraulis eurystole", # 1989, comes and goes, 5%, # ---- KEEP ---- ????
"Etropus cyclosquamus", # 1990, shelf flounder, jeff rester says bad # ---- TOSS ----
# "Etropus microstomus", # 1985, comes and goes, 2%, flounder, # ---- KEEP ----
# "Libinia dubia", # 1985, crab, 10%, # ---- KEEP ----
# "Libinia emarginata", # 1985, crab, 15%, # ---- KEEP ----
"Luidia clathrata", # 1990, sea star, 20%, # ---- TOSS ----
# "Mullus auratus", # 1985, goatfish, 10%, # ---- KEEP ----
# "Mustelus canis", # pres 1st yr, # ---- KEEP ----
# "Neobythites gilli", # 1985, brotula (fish), 5%, # ---- KEEP ----
# "Neomerinthe hemingwayi", # pres 1st yr, scorpionfish, 5%, # ---- KEEP ----
"Ogcocephalus declivirostris", # batfish rester says toss # ---- TOSS ----
# "Ophidion holbrookii", # 1999, eel, 5%, # ---- KEEP ----
# "Ovalipes floridaus", # 1985, crab, 5%, # ---- KEEP ----
# "Ovalipes stephensoni", # 1996, pres 1st yr, crab, # ---- KEEP ----
# "Paralichthys lethostigma", # southern flounder --- peculiar, # ---- KEEP ----
# "Parapenaeus politus", # 1987, shrimp, # ---- KEEP ----
"Pareques iwamotoi", # blackbar drump, rester says toss # ---- TOSS ----
# "Pareques umbrosus", # 1985, cubbyu (fish), 15%, # ---- KEEP ----
# "Parthenopoides massena", # 1987, a crab, 10%, # ---- KEEP ----
# "Peristedion gracile", # pres 1st yr, searobin ~10%, # ---- KEEP ----
"Pitar cordatus", # 1992, clam, 12%, # ---- TOSS ----
# "Pontinus longispinis", # pres 1st yr, scropionfish # ---- KEEP ----
"Reilla mulleri", # 1987, coral thing, 20%, # ---- TOSS ----
# "Rhomboplites aurorubens", # pres 1st yr, snapper, # ---- KEEP ----
# "Rimapenaeus constrictus", # 1988, 18%, shrimp # ---- KEEP ----
# "Rimapenaeus similis", # 1985, shrimp, 10%, # ---- KEEP ----
# "Sardinella aurita", # pres 1st year, sardinella, # ---- KEEP ----
# "Saurida caribbaea", # 1986, lizardfish, 5%, # ---- KEEP ----
# "Scomber japonicus", # pres 1st year, flickers out in 2000, chub mackerel, # ---- KEEP ----
# "Scomberomorus cavalla", # pres 1st year, 10%, king mackerel, # ---- KEEP ----
# "Seriola dumerili", # 1986, amberjack, # ---- KEEP ----
# "Sicyonia burkenroadi", # 1986, shrimp, # ---- KEEP ----
# "Soleocera vioscai", # 1988, shrimp, 10%, jumps to 40% in ~1995; # ---- KEEP ----
# "Squatina dumeril", # 1986, angelshark, # ---- KEEP ----
# "Squilla chydaea", # 1985, 15%, mantis shrimp, jump to ~80% 1990; hmm # ---- KEEP ---- ????
# "Squilla empusa", # 1985, 15%, mantis shrimp, jumps to 90% 1990; # ---- KEEP ---- ????
# "Steindachneria argentea", # pres 1st year, 4%, hake, # ---- KEEP ----
# "Symphurus civitatum", # 1985, 5%, tonguefish, comes and goes, # ---- KEEP ----
""
)

# ---- GoA ----
print(col_ext_full[reg=="goa", list(spp, orig_tax_id, full_firstYear, full_lastYear, full_propYear)], nrow=Inf)
man_remove_goa <- c(
"Actinauge verrilli", # 1999, 20%, anemone, high certainty; # ---- TOSS ----
"Alcyonidium pedunculatum", # 2007, bryozoan, 3%, low certainty; # ---- TOSS ----
# "Aphrocallistes vastus", # 1993 ('87 full), 5%, sponge, low certainty; # ---- KEEP ---- ????
"Aphrodita negligens", # 1999, 1%, worm thing, low certainty; # ---- TOSS ----
# "Arctomelon stearnsii", # 1987, 1%, snail, low certainty; # ---- KEEP ---- ????
"Asterias amurensis", # 1993, 5%, sea star, high certainty; # ---- TOSS ----
# "Bathyraja parmifera", # pres first yr, rare (1%), high certainty; # ---- KEEP ----
# "Brisaster latifrons", # 1993 ('87 full), 1%, sea star, high certainty; # ---- KEEP ----
"Buccinum plectrum", # 1990, 1%, whelk, low certainty; # ---- TOSS ----
"Buccinum scalariforme", #1996 ('93 full), 20%, whelk, high certainty; # ---- TOSS ----
"Ceramaster japonicus", # 1993, sea star, 5%, low certainty; # ---- TOSS ----
"Ceramaster patagonicus", # 1993, sea star, 15%, low certainty; # ---- TOSS ----
"Cheiraster dawsoni", # 1990, sea star, 1%, high certainty; # ---- TOSS ----
# "Chlamys rubida", # 1987, scallop, 1%, low certainty; # ---- KEEP ---- ????
# "Chorilia longipes", # 1984, 1%, crab, medium certainty; # ---- KEEP ----
"Chrysaora melanaster", # 1999, 1%, sea nettle, high certainty; # ---- TOSS ---- ????
"Clathrina blanca", # 1996, sponge, 5%, high certainty; # ---- TOSS ----
"Crossaster borealis", # 1993, sea star, 5%, high certainty; # ---- TOSS ----
"Crossaster papposus", # 1987, rare at firs tehn super common, high certainty # ---- TOSS ---- ????
# "Cucumaria fallax", # 1987, 1%, cucumber, high certainty; # ---- KEEP ---- ????
# "Cyanea capillata", # 1990 ('84 full), jellyfish, medium certainty; # ---- KEEP ----
"Diplopteraster multipes", # 1993 ('87 full), sea star, high certainty # ---- TOSS ----
"Dipsacaster borealis", # 1993, sea star, 10%, high certainty; # ---- TOSS ----
# "Doris odhneri", # 1999 ('87 full), slug, always rare (1-12%), high certainty; # ---- KEEP ---- ????
"Elassochirus cavimanus", # 1996, hermit, 10%, high certainty; # ---- TOSS ----
"Elassochirus gilli", # 1996, hermit, rare (1%), high certainty; # ---- TOSS ----
# "Enteroctopus dofleini", # 1987, giant octopus, 1%, high certainty; # ---- KEEP ---- ????
# "Eumicrotremus phrynoides", #1999, lumpsucker, 5%, medium certainty; # ---- KEEP ---- ????
"Evasterias troscheli", # 1993, sea star, 5%, low certainty; # ---- TOSS ----
"Fanellia compressa", # 1996, sea fan, 5%, low certainty; # ---- TOSS ----
# "Glebocarcinus oregonensis", # 1987, crab, 1%, high certainty; # ---- KEEP ---- ????
"Halichondria panicea", # 1993, sponge, 5%, high certainty; # ---- TOSS ----
"Halichondria sitiens", # 1996, sponge, 5%, low certainty; # ---- TOSS ----
"Halipteris willemoesi", # 2005, sea pen, low certainty; # ---- TOSS ----
# "Halocynthia aurantium", # 1987, sea peach, high certainty; # ---- KEEP ---- ????
"Halocynthia dumosa", # 1996, tunicate, 5%, low certainty; # ---- TOSS ----
"Halocynthia igaboja", # 2007, tunicate, 15%, low certainty; # ---- TOSS ----
"Henricia leviuscula", # 1993, sea star, 5%, low certainty; # ---- TOSS ----
# "Hippasteria phrygiana", # 1984 (only '87 in full, mltpl spp, updated?), 1% high certainty; # ---- KEEP ---- ????
"Hyas lyratus", # pres 1st yr, pres all yrs in full, 20%, high certainty; # ---- KEEP ----
# "Hydrolagus colliei", # flawed b/c of strat during col yr, high certainty # ---- KEEP ----
# "Labidochirus splendescens", # 1990, hermit, 1%, high certainty; # ---- KEEP ---- ????
"Laqueus californianus", # 1993 ('90 full), clam, low certainty; # ---- TOSS ----
# "Luidia foliolata", # 1987, sea star, 2%, high certainty; # ---- KEEP ---- ????
# "Lycodes diapterus", # 1984, ~5%, eelpout, high certainty; # ---- KEEP ----
# "Merluccius productus", # 1987, 5%, pacific hake, high certainty; # ---- KEEP ----
# "Mesocentrotus franciscanus", # 1987 ('84 full, but miss '87 full?), urchin, 1%, high certainty; # ---- KEEP ----
# "Metridium dianthus", # 1984, anemone, 2-20%, medium certainty; # ---- KEEP ----
"Metridium farcimen", # 1999, 1% then immediately 50%, anemone, high certainty; # ---- TOSS ----
# "Modiolus modiolus", # 1987, mussel, 1%-20%, high certainty; # ---- KEEP ---- ????
"Mycale loveni", # 1993, tree sponge, 15%, low certainty; # ---- TOSS ----
# "Myoxocephalus jaok", # 1987, sculpin, 3%, high certainty; # ---- KEEP ---- ????
"Myxilla incrustans", # 1990, sponge, 3%, low certainty; # ---- TOSS ----
"Neoesperiopsis infundibula", # 1996, sponge, 2%, low certaint # ---- TOSS ----
"Neptunea amianta", # 1993, whelk, 2%, low certainty, # ---- TOSS ----
# "Neptunea lyrata", # 1990 ('87 full), 5%, low certainty; # ---- KEEP ---- ????
# "Neptunea pribiloffensis", # 1987 ('84 full), 2%, low certainty; # ---- KEEP ----
"Ophiopholis aculeata", # 1993, star, 2%, high certainty; # ---- TOSS ----
"Ophiura sarsii", # 1990, star, 1%, high certainty; # ---- TOSS ----
"Orthasterias koehleri", # 1993, star, 10%, high certainty; # ---- TOSS ----
# "Pagurus aleuticus", # 1990 ('87 full), 5%, hermit, high certainty; # ---- KEEP ---- ????
# "Pagurus brandti", # 1990 ('87 full), 1%, hermit, low certainty; # ---- KEEP ---- ????
"Pagurus capillatus", # 1996 ('90 full), 5%, hermit, low certainty; # ---- TOSS ----
"Pagurus confragosus", # 1993, hermit, 1%, medium certainty; # ---- TOSS ----
"Pagurus kennerlyi", # 1996, 5%, hermit, low certainty; # ---- TOSS ----
# "Pagurus ochotensis", # 1987, hermit, 2%, medium certainty; # ---- KEEP ---- ????
"Pagurus trigonocheirus", # 1990, 2%, low certainty; # ---- TOSS ----
# "Pandalus jordani", # 1987, 15%, shrimp, low certainty; # ---- KEEP ---- ????
# "Pandalus tridens", # 1990, 1%, shrimp, medium certainty; # ---- KEEP ---- ????
"Phacellophora camtschatica", # 1999, 1%, jelly, high certainty; # ---- TOSS ---- ????
"Pododesmus macrochisma", # 1990, 1%, jingle, high certainty; # ---- TOSS ---- ????
"Pseudarchaster parelii", # 1993, sea star, 20%, low certainty; # ---- TOSS ----
"Pseudostichopus mollis", # 1999 ('96 full), cucumber, 5%, high certainty; # ---- TOSS ----
"Pteraster militaris", # 1993, sea star, 5%, low certainty; # ---- TOSS ----
"Pteraster tesselatus", # 1987, sea star, low certainty, 1-20%, # ---- TOSS ----
"Ptilosarcus gurneyi", # 1993, sea pen, high certainty, 1%, # ---- TOSS ----
"Pyrulofusus harpa", # 1990, whelk, 1%-5%, high certainty; # ---- TOSS ---- ????
# "Rhinolithodes wosnessenskii", # miss 1987, but pres all yrs full; crab; 1%, high cert; # ---- KEEP ----
# "Sebastes elongatus", # only miss in yr where eastern part not sampled # ---- KEEP ----
# "Sebastes helvomaculatus", # same situation as above (miss b/c of east not sampled) # ---- KEEP ----
"Sebastes melanostictus", # was grouped with rougheye rockfish originally; # ---- TOSS ----
"Solaster dawsoni", # 1993 sea star # ---- TOSS ----
"Solaster endeca", # 1993 ('87 full), sea star, 5%, low certainty; # ---- TOSS ---- ????
"Stegophiura ponderosa", # 1993 sea star; 1% then quickly 25%; high certainty; # ---- TOSS ----
# "Stenobrachius leucopsarus", # 1990 (pres 1st yr); lampfish, 4%, medium certainty; # ---- KEEP ----
# "Strongylocentrotus fragilis", # pres 1st yr, urchin, medium certainty; # ---- KEEP ----
# "Strongylocentrotus pallidus", # 1987, urchin, 1%, comes and goes lots; low certainty; # ---- KEEP ---- ????
"Strongylocentrotus polyacanthus", # 1999, urchin, lots come and go, but prob not ID'd early; # ---- TOSS ---- ????
"Styela rustica", # 1993, sea potato, 1%, high certainty; # ---- TOSS ----
# "Stylasterias forreri", # pres first year, high certainty, crab, ~5%; # ---- KEEP ----
"Suberites domuncula", # 2007, sponge, low certainty, 10%, # ---- TOSS ----
"Suberites ficus", # 1990 ('87 full), sponge, 5%, low certainty; # ---- TOSS ---- ????
"Synallactes challengeri", # 2001, worm thing, 5%, high certainty; # ---- TOSS ----
"Terebratalia transversa", # 1987, brachiopod, 1%, low certainty; # ---- TOSS ---- ????
""
)


# ---- NEUS ----
# print(col_ext_full[reg=="neus", list(spp, orig_tax_id, full_firstYear, full_lastYear, full_propYear)], nrow=Inf)
man_remove_neus <- c(
# "Acipenser oxyrinchus",
# "Anarhichas lupus",
# "Anchoa hepsetus",
# "Antigonia capros",
# "Argentina silus",
# "Argentina striata",
# "Atlantopandalus propinqvus",
# "Bairdiella chrysoura",
# "Bathynectes longispina",
# "Bathypolypus arcticus",
# "Brevoortia tyrannus",
# "Carcharhinus plumbeus",
# "Centropristis philadelphica",
# "Chaceon quinquedens",
# "Chionoecetes opilio",
# "Chlorophthalmus agassizi",
# "Conger oceanicus",
# "Crangon septemspinosa",
# "Cryptacanthodes maculatus",
# "Cyclopterus lumpus",
# "Cynoscion regalis",
# "Dasyatis centroura",
# "Dasyatis say",
# "Decapterus punctatus",
# "Dichelopandalus leptocerus",
# "Dipturus laevis",
# "Etropus microstomus",
# "Etrumeus teres(sadina)",
# "Foetorepus agassizii",
# "Gymnura altavela",
# "Lagodon rhomboides",
# "Larimus fasciatus",
# "Lebbeus polaris",
# "Leiostomus xanthurus",
# "Leptoclinus maculatus",
# "Liparis atlanticus",
# "Lithodes maja",
# "Lolliguncula brevis",
# "Lopholatilus chamaeleonticeps",
# "Lumpenus lampretaeformis",
# "Lycenchelys verrillii",
# "Macrorhamphosus scolopax",
# "Maurolicus weitzmani",
# "Melanostigma atlanticum",
# "Menticirrhus americanus",
# "Menticirrhus saxatilis",
# "Micropogonias undulatus",
# "Monolene sessilicauda",
# "Morone saxatilis",
# "Mustelus canis",
# "Myliobatis freminvillii",
# "Myoxocephalus aenaeus",
# "Myxine glutinosa",
# "Nemichthys scolopaceus",
# "Ophidion marginatum",
# "Osmerus mordax",
# "Ovalipes ocellatus",
# "Ovalipes stephensoni",
# "Pandalus borealis",
# "Pandalus montagui",
# "Parasudis truculenta",
# "Pasiphaea multidentata",
# "Peristedion miniatum",
# "Petromyzon marinus",
# "Phycis chesteri",
# "Polymixia lowei",
# "Pomatomus saltatrix",
# "Pontinus longispinis",
# "Pontophilus norvegicus",
# "Prionotus alatus",
# "Scyliorhinus retifer",
# "Sicyonia brevirostris",
# "Syacium papillosum",
# "Symphurus plagiusa",
# "Synagrops bellus",
# "Syngnathus fuscus",
# "Synodus foetens",
# "Synodus poeyi",
# "Trachinocephalus myops",
# "Trachurus lathami",
# "Trichiurus lepturus",
# "Triglops murrayi",
# "Trinectes maculatus",
# "Ulvaria subbifurcata",
# "Zenopsis conchifera",
""
)


# ---- NEWF ----
print(col_ext_full[reg=="newf", list(spp, orig_tax_id, full_firstYear, full_lastYear, full_propYear)], nrow=Inf)
man_remove_newf <- c(
# "Amblyraja jenseni",
# "Arctozenus risso kroyeri",
# "Artediellus atlanticus",
# "Ceramaster granularis",
# "Cross papposus",
# "Eumicrotremus spinosus",
# "Flabellum macandrewi",
# "Gersemia rubiformis",
# "Hyas araneus",
# "Liparis atlanticus",
# "Lycodes reticulatus",
# "Myoxocephalus scorpioides",
# "Myxine glutinosa",
# "Notacanthus chemnitzii",
# "Ophiura robusta",
# "Paragorgia arborea",
# "Pasiphaea princeps",
# "Pasiphaea tarda",
# "Polyacanthonotus rissoanus",
# "Poromitra capito",
# "Primnoa resedaeformis",
# "Sergestes arcticus",
# "Triglops murrayi",
""
)


# ---- SA ----
# print(col_ext_full[reg=="sa", list(spp, orig_tax_id, full_firstYear, full_lastYear, full_propYear)], nrow=Inf)
man_remove_sa <- c(
# "Acanthostracion quadricornis",
# "Alectis ciliaris",
# "Aluterus schoepfii",
# "Ariopsis felis",
# "Bagre marinus",
# "Calamus leucosteus",
# "Calappa flammea",
# "Carcharhinus acronotus",
# "Carcharhinus limbatus",
# "Centropristis philadelphica",
# "Citharichthys spilopterus",
# "Diplectrum formosum",
# "Dyspanopeus sayi",
# "Etrumeus teres(sadina)",
# "Farfantepenaeus duorarum",
# "Haemulon aurolineatum",
# "Harengula jaguana",
# "Hippocampus erectus",
# "Hypleurochilus geminatus",
# "Menticirrhus saxatilis",
# "Mobula hypostoma",
# "Mustelus canis",
# "Myliobatis freminvillii",
# "Paralichthys squamilentus",
# "Pilumnus sayi",
# "Prionotus rubio",
# "Rachycentron canadum",
# "Raja eglanteria",
# "Rhinobatos lentiginosus",
# "Rhinoptera bonasus",
# "Sardinella aurita",
# "Sphyrna lewini",
# "Syacium papillosum",
# "Symphurus plagiusa",
# "Xiphopenaeus kroyeri",
""
)


# ---- Shelf ----
print(col_ext_full[reg=="shelf", list(spp, orig_tax_id, full_firstYear, full_lastYear, full_propYear)], nrow=Inf)
man_remove_shelf <- c(
# "Alosa pseudoharengus",
# "Alosa sapidissima",
# "Arctozenus risso",
# "Artediellus atlanticus",
# "Artediellus uncinatus",
# "Aspidophoroides monopterygius",
# "Citharichthys arctifrons",
# "Cyclopterus lumpus",
# "Enchelyopus cimbrius",
# "Eumicrotremus spinosus",
# "Helicolenus dactylopterus",
# "Leptagonus decagonus",
# "Leptoclinus maculatus",
# "Leucoraja erinacea",
# "Lumpenus lampretaeformis",
# "Lycodes vahlii",
# "Mallotus villosus",
# "Myxine glutinosa",
# "Peprilus triacanthus",
# "Reinhardtius hippoglossoides",
# "Urophycis chuss",
""
)


# ---- WCAnn ----
print(col_ext_full[reg=="wcann", list(spp, orig_tax_id, full_firstYear, full_lastYear, full_propYear)], nrow=Inf)
man_remove_wcann <- c(
# "Anthopleura xanthogrammica",
# "Bathylagus pacificus",
# "Cataetyx rubrirostris",
# "Chesnonia verrucosa",
# "Chrysaora fuscescens",
# "Citharichthys xanthostigma",
# "Dosidicus gigas",
# "Halocynthia igaboja",
# "Heteropolypus ritteri",
# "Metacarcinus gracilis",
# "Mustelus henlei",
# "Neognathophausia ingens",
# "Neverita lewisii",
# "Octopus rubescens",
# "Philine bakeri",
# "Pyrosoma atlanticum",
# "Sergestes similis",
# "Stenobrachius leucopsarus",
# "Theragra chalcogramma",
# "Thetys vagina",
# "Trachurus symmetricus",
# "Urticina columbiana",
""
)


# ---- WCTri ----
print(col_ext_full[reg=="wctri", list(spp, orig_tax_id, full_firstYear, full_lastYear, full_propYear)], nrow=Inf)
man_remove_wctri <- c(
"Apostichopus leukothele", # 1998, 30%, cucumber; # ---- TOSS ----
# "Argentina sialis", # 1983, 5%, argentine, # ---- KEEP ---- ????
# "Bathyagonus pentacanthus", # 1983, poacher, 5%, # ---- KEEP ----
"Brisaster latifrons", # 1989, urchin, 5%, then BOOM 60%; # ---- TOSS ----
# "Cymatogaster aggregata", # 1983, shiner perch, 5%, # ---- KEEP ----
# "Enteroctopus dofleini", # 1989, 10%, ocotopus; # ---- KEEP ----
# "Hexagrammos decagrammus", # 1980, greenling, 2%, # ---- KEEP ----
# "Hippasteria phrygiana", # sea star, 5%, 1995 # ---- KEEP ---- ????
"Lepidopsetta bilineata", # rock sole --- probably previous confused with polyxstra like in alaska # ---- TOSS ----
"Liponema brevicorne", # 1998, anemone 10% # ---- TOSS ----
# "Luidia foliolata", # 1983, sand star, comes and goes, 20% - 80%+ # ---- KEEP ---- ????
# "Lycodes diapterus", # 1980, eelpout, 2%, # ---- KEEP ----
# "Lycodes pacificus", # 1980, eelpout, 5%, # ---- KEEP ----
# "Mediaster aequalis", # 1995, sea star, 5%, # ---- KEEP ---- ????
# "Metridium dianthus", # 1980, anemone, 1%, # ---- KEEP ----
"Metridium farcimen", # 2001, giant anemone, 40%, # ---- TOSS ----
# "Pandalus eous", # present 1st yr, shrimp, # ---- KEEP ----
# "Peprilus simillimus", # 1983, pompano, # ---- KEEP ----
# "Pisaster brevispinus", # 1983, comes and goes, 2%, # ---- KEEP ----
# "Platichthys stellatus", # 1980, starry flounder # ---- KEEP ----
# "Pleurobranchaea californica", # 1989, nudibranch 8% # ---- KEEP ----
# "Pleuronichthys verticalis", # pres 1st yr, flatfish, 2%, # ---- KEEP ----
# "Psettichthys melanostictus", # 1980, 5%, flatfish # ---- KEEP ----
# "Pycnopodia helianthoides", # 1986, 5% for a while, them booms to 50% # ---- KEEP ---- ????
# "Raja inornata", # 1980, ray, 5%, # ---- KEEP ----
# "Rathbunaster californicus", # 1995, sea star, 15%, goes to 40% # ---- KEEP ----
# "Rossia pacifica", # pres 1st yr, octopus, 5%, # ---- KEEP ----
"Sardinops sagax", # 1992, pichard, 30%, # ---- TOSS ----
# "Scomber japonicus", # 1983, mackerel, 2%, # ---- KEEP ----
# "Sebastes brevispinis", # leaver, # ---- KEEP ----
# "Sebastes semicinctus", # 1989, 10%, rockfish, # ---- KEEP ---- ????
# "Strongylocentrotus fragilis", # 1986, urchin, 10%, # ---- KEEP ---- ????
"Stylasterias forreri", # 1995, sea star, 10% # ---- TOSS ----
"Tritonia diomedea", # 1989, slug, 10%, comes and goes a bit, but booms later # ---- TOSS ---- ????
# "Xeneretmus latifrons", # 1983, poacher, 5%, # ---- KEEP ----
# "Zaniolepis latipinnis", # 1983, combfish, 8%, # ---- KEEP ----
""
)


man_remove <- list(
	ebs = man_remove_ebs,
	ai = man_remove_ai,
	goa = man_remove_goa,
	wctri = man_remove_wctri,
	wcann = man_remove_wcann,
	gmex = man_remove_gmex,
	sa = man_remove_sa,
	neus = man_remove_neus,
	shelf = man_remove_shelf,
	newf = man_remove_newf
)

save(man_remove, file="~/Documents/School&Work/pinskyPost/trawl/trawlDiversity/pkgBuild/spp_check/man_remove.RData")