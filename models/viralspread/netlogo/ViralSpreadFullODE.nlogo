;; ----- DECLARATIONS -----
;; use arrays
extensions [array]

;; Create "breed" of 'turtle' called cell
breed [cells cell]          ;; eventually separate mesophyll?
breed [vasculars vascular]  ;; phloem vascular bundles 

;; Create "breed" of 'link' called phloem
directed-link-breed [phloems phloem]  
;; Note: I know there shouldn't be an s, but it wanted a different plural :(
undirected-link-breed [plasmodesmata plasmodesma]

;; Declare global variables
globals
[
  num-infected   ;; keep track of number infected
  num-viruses    ;; keep track of the viral particles
  mod-num-viruses  ;; keep track of the modified viral particles
  lysed-cells    ;; keep track of apoptotic cells
  infect-per-leaf ;; array of infection counts
]


;; Declare the cell-specific (turtle breed) variables
turtles-own
[
  infected?         ;; if true, the cell is infectious
  resistant?        ;; if true, cell can't be infected
  num-plasmodesmata ;; number of connections with other cells
  leaf              ;; for IC setup
  
  dna-gap           ;; track the number of genomes in the nucleus
  dna-ccc           ;;
  mod-dna-gap       ;;
  mod-dna-ccc       ;;
  rna-19s           ;;
  rna-35s           ;;
  mod-rna-35s       ;;
  protein-one       ;;
  protein-two       ;;
  protein-three     ;;
  protein-four      ;;
  protein-four-sub  ;;
  protein-five      ;;
  protein-six       ;;
  int-virions       ;;
  viral-count       ;;
  mod-int-virions   ;;
  mod-viral-count   ;;
  sar-level         ;; amount of "salicylic acid" present in the cell
]

;; Setup procedures
to setup
  clear-all    ;; remove anything from previous runs
  setup-stem   ;; make a stem structure to connect the leaves
  setup-cells  ;; set up the cells
  setup-leaf   ;; set up the connections between cells
  set infect-per-leaf array:from-list n-values num-leaves [0] ;; initialize as 1xn array of 0's
  ask n-of initial-infection-sites cells with [leaf = 1]  ;; infect initial cells
                                         ;; currently all in leaf 1, but could change
    [ 
      become-infected 
      set dna-gap founder-population-viruses ;; start with x genomes per infected cell
    ]
  ask links [ set color white ]  ;; make the symplastic connections white
  ask phloems 
    [ 
      set color green            ;; make the vascular bundles green and thicker
      set thickness 0.3
    ]
  reset-ticks  ;; reset timer from previous run
  if file-exists? "netlogo_sim.csv" [file-delete "netlogo_sim.csv"]  ;; clean slate for output
end

to setup-stem
  set-default-shape vasculars "plant" ;; uses stem-looking shape for vasculature
  ;; Create and position vascular cells
  create-vasculars 1 [ setxy 0 -16 ]
  create-vasculars 1 
  [ 
    setxy 0 -6 
    create-phloem-from vascular 0 ;; directed stem connection
  ]
  create-vasculars 1 
  [ 
    setxy 0 6 
    create-phloem-from vascular 1 ;; directed stem connection
  ]
  create-vasculars 1 
  [ 
    setxy 0 16 
    create-phloem-from vascular 2 ;; directed stem connection
  ]
  ask vasculars         ;; make the vasculature stand out from the leaf cells
    [ 
      set color green 
      set size 2.5
    ]
end

to setup-cells
  set-default-shape cells "circle"   ;; uses circle shape for displayed cells
  let cells-per-leaf floor (num-cells / num-leaves)
  let extra-cells remainder num-cells num-leaves
  
  let angular-spacing 10 ;; degrees between leaves
  if angular-spacing * num-leaves > 360   ;; throw an error if the leaves don't fit
   [ error "Not enough space for leaves"]
  let leaf-angle (360 - angular-spacing * num-leaves) / num-leaves ;; angular width per leaf
  let rmin 0.1 * max-pxcor  ;; minimum distrance from the origin
  let rmax min list max-pxcor max-pycor  ;; max distance -- 
  
  foreach n-values  extra-cells [?]  ;; ? here is just identity, making n-values 0:n-1
  [create-cells cells-per-leaf + 1 ;; 
    [
    set size 0.5            ;; netlogo is a freak, so '?' == loop iteration - 1
    let theta (? * (leaf-angle + angular-spacing)) + random-float leaf-angle  
    let r random-radius rmin rmax
    let x r * cos(theta)
    let y r * sin(theta)
    setxy x y
    set leaf ? + 1
    become-susceptible
    ]
   ]   
  ;; now do it again for leaves without extra cells 
  foreach n-values  (num-leaves - extra-cells) [? + extra-cells]  ;; ? here is just identity, making n-values 0:n-1
  [create-cells cells-per-leaf ;; 
    [
    set size 0.5            ;; netlogo is a freak, so '?' == loop iteration - 1
    let theta (? * (leaf-angle + angular-spacing)) + random-float leaf-angle  
    let r random-radius rmin rmax
    let x r * cos(theta)
    let y r * sin(theta)
    setxy x y
    set leaf ? + 1
    become-susceptible
    ]
   ] 
end

to setup-leaf
  ask cells 
  [ 
    ;; using random-possion to give a more realistic distribution of links per cell
    let my-leaf leaf
    while [ count my-links < random-poisson avg-num-plasmodesmata ]
    [ 
      let choice (min-one-of (other cells with [(not link-neighbor? myself) and leaf = my-leaf])
                              [distance myself])
      if choice != nobody [ create-plasmodesma-with choice ]
    ]
  ]
  ;; to make the network look nicer
;;repeat 10 [ layout-spring turtles links 0.3 (world-width / (sqrt num-cells)) 1 ]
  ;; connect leaves to the vasculature
  ask vascular 0 
  [ 
    let choice ( one-of (cells with [leaf = 1]) )
    if choice != nobody [ create-phloem-from choice ]
    become-susceptible
  ]
  ask vascular 1 
  [ 
    let choice ( one-of (cells with [pxcor > 0 and pycor < 0]) )
    if choice != nobody [ create-phloem-from choice ]
    become-susceptible
  ]
  ask vascular 2 
  [ 
    let choice ( one-of (cells with [pxcor < 0 and pycor > 0]) )
    if choice != nobody [ create-phloem-from choice ]
    become-susceptible
  ]
  ask vascular 3
  [ 
    let choice ( one-of (cells with [pxcor > 0 and pycor > 0]) )
    if choice != nobody [ create-phloem-from choice ]
    become-susceptible
  ]
end


to-report random-radius [rmin rmax]
  ;; picks a radius at random such that points will be uniform in a circle
  ;; pdf is r/const ; use inverse transform method
  let u random-float 1
  let r sqrt(u * (rmax ^ 2 - rmin ^ 2) + rmin ^ 2 )
  report r
end


;; ----- MAIN -----

;; This is the main process - what is activated when the user clicks "go"
to go
  ;; if the user specified a duration, stop after that
  if specified-duration
  [
    if ticks > num-ticks [ stop ] 
  ]
  ;; otherwise, check if all cells are infected
  if all? cells [ infected? ] 
  [  
    ;; if they are, run 40 more times (make the graph look nice)
    let i 0 
    while [i <= 40] 
    [
      spread-virus
      set i i + 1
      tick
    ]
    stop  ;; stop model completely
  ]
  
  ;; otherwise, continue to assemble and spread the virus
  assemble-virus
  spread-virus
  ;;spread-sar
  ;;do-apoptosis-checks
  record-data
  tick
end

;; ----- DATA OUTPUT -----

to record-data
  file-open "netlogo_sim.csv"
  let leaf-infect-list array:to-list infect-per-leaf
  foreach leaf-infect-list  ;; iterate through this list
  [ file-write ?            ;; type number then comma
    file-type ","]
  file-type "\n"            ;; newline
  file-close                ;; apparently bad things happen if you don't
end


;;  ----- CELL PROCEDURES -----

;; Cell Procedures
to become-infected
  set infected? true
  set resistant? false
  set color red
  ;; determine number of viruses entering the susceptible cell
  ;let new-viruses ( sum [viral-count] of link-neighbors )
  set dna-gap dna-gap + 1
  set num-viruses num-viruses + 1
  array:set infect-per-leaf (leaf - 1) (array:item infect-per-leaf (leaf - 1) + 1) ;; increment counter
end

to mod-become-infected
  set infected? true
  set resistant? false
  set color red
  set mod-dna-gap mod-dna-gap + 1
  set mod-num-viruses mod-num-viruses + 1
end

to become-susceptible
  set infected? false
  set resistant? false
  set color green
  set dna-gap 0          ;; track the number of genomes in the nucleus
  set dna-ccc 0          ;;
  set mod-dna-gap 0      ;;
  set mod-dna-ccc 0      ;;
  set rna-19s 0          ;;
  set rna-35s 0          ;;
  set mod-rna-35s 0      ;;
  set protein-one 0      ;;
  set protein-two 0      ;;
  set protein-three 0    ;;
  set protein-four 0     ;;
  set protein-four-sub 0 ;;
  set protein-five 0     ;;
  set protein-six 0      ;;
  set int-virions 0      ;;
  set viral-count 0      ;;
  set mod-int-virions 0  ;;
  set mod-viral-count 0  ;;
end

;; ---- SAR Molecule Spread & Apoptosis ----

;; Increase levels of signalling molecule based on neighbours' levels
;; Assuming the infected cells are unable to produce sar signal molecules, but 
to spread-sar
  ;; initial generation of sar based on neighbouring cells being infected
  ask cells
  [
    let neighboured false           ;; remember if a neighbour is infected
    ask link-neighbors 
    [
      if (infected?) [ set neighboured true ]
    ]
    ;; if a neighbour is infected, increase this cell's sar signalling molecule
    if (neighboured) [ set sar-level sar-level + 1 ]
  ]
  ;; assume the resistant cells are better able to spread the sar molecule
  ask cells with [resistant?]
  [
    let shared-sar sar-level / 2     ;; saying half the sar-level will be spread
    ask link-neighbors 
    [ 
      set sar-level sar-level + shared-sar 
      if sar-level >= resistance-threshold [ set resistant? true ]
    ]
    set sar-level sar-level + shared-sar / 2  ;; not much kept in the cell
  ]
  ;; assume the susceptible cells can produce an okay amount of sar molecule
  ask cells with [not infected? and not resistant?]
  [
    let shared-sar sar-level / 5     ;; fifth of the molecules will be shared
    ask link-neighbors
    [
      set sar-level sar-level + shared-sar
      if sar-level >= resistance-threshold [ set resistant? true ]
    ]
    set sar-level sar-level + shared-sar / 2
  ]
end  

;; Determine whether the cell can destroy itself to help prevent viral spread
to do-apoptosis-checks
  ;; if the cell is infected and can still lyse, check 
  ask cells with [infected?] 
  [
    if (sar-level >= lysis-threshold)
      [ 
        set color violet
        die 
        set lysed-cells lysed-cells + 1
      ] 
  ]
end


;;  ----- VIRUS PROCEDURES -----

;; Procedure governing spread to neighbouring cells
to spread-virus
  ask turtles with [infected?]
    [ 
      ;;if (viral-count + mod-viral-count > 5)
        ;;[
          ask link-neighbors with [not resistant? and not infected?]
            [ 
              if (random-float 100 < viral-spread-chance)
                [ 
                  become-infected                   ;; If chance has it, infect the cell
                  set num-infected num-infected + 1 ;; increase the count for infected cells
                ] 
            ] 
        ;;]
    ]
end

to assemble-virus
  
  ;; Genome parameters
  let max-genomes-in-nucleus 100
  let intracellular-reinfection-rate 0.01
  let gapped-dna-repair-rate 0.1
  let dna-modification-rate 0.01
  let dna-degradation-rate 0.00001
  
  ;; RNA parameters
  let transcription-rate-19s 0.05
  let degradation-rate-19s 0.001155
  let transcription-rate-35s 0.0653
  let degradation-rate-35s 0.001155
  let packaging-rate 0.1
  let frac-unspliced 0.3
  
  ;; Protein parameters
  let translation-rate 0.1
  let degradation-rate 0.0001
  let p4-splicing-rate 1
  
  ;; Virion parameters
  let anchoring-rate 1
  let virion-degradation-rate 0.01
  let virion-exit-rate 0.1
  
  ;; Time step
  let delta-t 0.001
  
  ask cells with [infected?]
    [
      set dna-gap          dna-gap + delta-t * ( intracellular-reinfection-rate * viral-count * ( max-genomes-in-nucleus - dna-gap - mod-dna-gap ) - gapped-dna-repair-rate * dna-gap - dna-modification-rate * dna-gap )
      set dna-ccc          dna-ccc + delta-t * ( gapped-dna-repair-rate * dna-gap - dna-degradation-rate * dna-ccc - dna-modification-rate * dna-ccc )
      set mod-dna-gap      mod-dna-gap + delta-t * ( intracellular-reinfection-rate * mod-viral-count * ( max-genomes-in-nucleus - dna-gap - mod-dna-gap ) - gapped-dna-repair-rate * mod-dna-gap + dna-modification-rate * dna-gap )
      set mod-dna-ccc      mod-dna-ccc + delta-t * ( gapped-dna-repair-rate * mod-dna-gap - dna-degradation-rate * mod-dna-ccc + dna-modification-rate * dna-ccc )
      set rna-19s          rna-19s + delta-t * ( transcription-rate-19s * dna-ccc - degradation-rate-19s * rna-19s )
      set rna-35s          rna-35s + delta-t * ( transcription-rate-35s * dna-ccc - degradation-rate-35s * rna-35s - packaging-rate * protein-four-sub * protein-five * rna-35s * frac-unspliced )
      set mod-rna-35s      mod-rna-35s + delta-t * ( transcription-rate-35s * mod-dna-ccc - degradation-rate-35s * mod-rna-35s - packaging-rate * protein-four-sub * protein-five * mod-rna-35s * frac-unspliced )
      set protein-one      protein-one + delta-t * ( translation-rate * ( rna-35s + mod-rna-35s ) * frac-unspliced - degradation-rate * protein-one )
      set protein-two      protein-two + delta-t * ( translation-rate * ( rna-35s + mod-rna-35s ) * frac-unspliced - degradation-rate * protein-two )
      set protein-three    protein-three + delta-t * ( translation-rate * ( rna-35s + mod-rna-35s ) - degradation-rate * protein-three - anchoring-rate * protein-three * ( int-virions + mod-int-virions ) )
      set protein-four     protein-four + delta-t * ( translation-rate * ( rna-35s + mod-rna-35s ) - degradation-rate * protein-four - p4-splicing-rate * protein-four )
      set protein-four-sub protein-four-sub + delta-t * ( p4-splicing-rate * protein-four - degradation-rate * protein-four-sub - packaging-rate * protein-four-sub * protein-five * ( rna-35s + mod-rna-35s ) * frac-unspliced )
      set protein-five     protein-five + delta-t * ( translation-rate * ( rna-35s + mod-rna-35s ) - degradation-rate * protein-five - packaging-rate * protein-four-sub * protein-five * ( rna-35s + mod-rna-35s ) * frac-unspliced )
      set protein-six      protein-six + delta-t * ( translation-rate * rna-19s - degradation-rate * protein-six )
      set int-virions      int-virions + delta-t * ( packaging-rate * protein-four-sub * protein-five * rna-35s * frac-unspliced - anchoring-rate * protein-three * int-virions )
      set viral-count      viral-count + delta-t * ( anchoring-rate * protein-three * int-virions - intracellular-reinfection-rate * viral-count - virion-degradation-rate * viral-count - virion-exit-rate * viral-count )
      set mod-int-virions  mod-int-virions + delta-t * ( packaging-rate * protein-four-sub * protein-five * mod-rna-35s * frac-unspliced - anchoring-rate * protein-three * mod-int-virions )
      set mod-viral-count  mod-viral-count + delta-t * ( anchoring-rate * protein-three * mod-int-virions - intracellular-reinfection-rate * mod-viral-count - virion-degradation-rate * mod-viral-count - virion-exit-rate * mod-viral-count )
    ]
    
    ;; oh thank god it's over
    
end
@#$#@#$#@
GRAPHICS-WINDOW
238
10
672
465
16
16
12.85
1
10
1
1
1
0
0
0
1
-16
16
-16
16
0
0
1
ticks
30.0

BUTTON
116
288
179
321
go
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
22
102
196
135
initial-infection-sites
initial-infection-sites
1
4
4
1
1
NIL
HORIZONTAL

SLIDER
20
10
192
43
num-cells
num-cells
1
1000
217
1
1
NIL
HORIZONTAL

SLIDER
20
49
208
82
avg-num-plasmodesmata
avg-num-plasmodesmata
1
10
4
1
1
NIL
HORIZONTAL

PLOT
231
475
454
672
#infected over time
ticks
#infected cells
0.0
50.0
0.0
50.0
true
false
"set-plot-x-range 0 num-ticks\nset-plot-y-range 0 num-cells" ""
PENS
"infected" 1.0 0 -13791810 true "" "plot num-infected"

BUTTON
38
289
101
322
setup
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
25
356
197
389
viral-spread-chance
viral-spread-chance
1
20
15
0.5
1
NIL
HORIZONTAL

SLIDER
21
227
193
260
num-ticks
num-ticks
0
500
102
1
1
NIL
HORIZONTAL

SWITCH
22
190
177
223
specified-duration
specified-duration
0
1
-1000

PLOT
466
485
666
635
Number of Viruses
ticks
#viruses
0.0
10.0
1.0
1000.0
true
false
"set-plot-x-range 0 num-ticks\nset-plot-y-range 0 num-cells" ""
PENS
"default" 1.0 0 -16777216 true "" "plot (num-viruses + mod-num-viruses)"

SLIDER
22
141
215
174
founder-population-viruses
founder-population-viruses
1
30
10
1
1
NIL
HORIZONTAL

SLIDER
24
446
196
479
lysis-threshold
lysis-threshold
1
10000
10000
1
1
NIL
HORIZONTAL

PLOT
19
490
219
640
Lysed Cells
ticks
#cells
0.0
10.0
0.0
10.0
true
false
"set-plot-x-range 0 num-ticks\nset-plot-y-range 0 10" ""
PENS
"default" 1.0 0 -16777216 true "" "plot lysed-cells"

SLIDER
23
410
195
443
resistance-threshold
resistance-threshold
1
100
94
1
1
NIL
HORIZONTAL

SLIDER
206
314
378
347
num-leaves
num-leaves
1
10
3
1
1
NIL
HORIZONTAL

@#$#@#$#@
## WHAT IS IT?

This is a model tracking viral spread in a plant-structure made of leaf cells and vascular cells. The spatial locations of the cells are divided into "leaves", and individual cells are connected by plasmodesmata. 


## HOW IT WORKS

The model begins by creating all of the cells - it places them in space (in 4 leaves separated by the coodrinate axes) and iteratively goes through the cells, giving each a random number of connections. 

All the cells start as susceptible, except for some small initial infection. Then the virus spreads along the plasmodesmata connections as time progresses.

In this model there's also a chance of cells becoming resistant. This is triggered by either proximity to infected cells, or by spreading among 


## HOW TO USE IT

Set the sliders as you see fit and click setup. Hit go to see what happens!


## THINGS TO NOTICE

How does the viral infection curve change as # of connections is increased and viral spread chance is decreased?


## THINGS TO TRY

Move the sliders! And note the variations in the random setup, especially as it concerns connections between different leaves.


## EXTENDING THE MODEL

Giving appropriate parameters for the viral spread chance. 


## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)


## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)


## CREDITS AND REFERENCES

UWaterloo's 2015 iGEM team 
Zoë Humphries & Robert Gooding-Townsend
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270

@#$#@#$#@
NetLogo 5.2.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180

@#$#@#$#@
0
@#$#@#$#@
