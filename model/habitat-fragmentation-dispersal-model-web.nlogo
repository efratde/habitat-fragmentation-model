;; =====================================================================================
;; Habitat-fragmentation / relatedness-dependent dispersal model — NetLogo Web (browser) DEMO
;; REDUCED variant of the desktop model `habitat-fragmentation-dispersal-model.nlogo`.
;; Same simulation engine; the `rnd` weighted sampling is reimplemented in pure NetLogo
;; (see `wn-with-repeats`). The `stats`/`profiler`/`rnd`/`array` extensions are NOT used here
;; (unsupported on NetLogo Web), so the 8 stats-based regression/quantile analysis plots,
;; BehaviorSpace experiments and the profiler button are omitted. See the Info tab for details
;; and use the desktop model for the full statistical-analysis version.
;; =====================================================================================
extensions [table csv]
breed [ plants plant ] ;; plants haing a relatedness-dependent dispersal strategy. Can only reproduce with other rd.plants
breed [ seeds seed ]         ;; seeds of both types


breed [sides side]     ; the four sides of the selection rectangle


globals [ fixation_p htz_dispersal patch_number mutations.overall fixation.map pollination.distances multi_year_data multi_year_data_neg multi_year_data_neu multi_year_data_pos
  current-state ; "not-started", "selecting", "waiting-to-drag", "dragging"
  select-x      ; coordinates for the start of the select box
  select-y
  drag-x        ; coordinates for the start of a drag operation
  drag-y
  selected      ; agentset of currently selected circles
  plants.all plants.all.0 plants.all.1 plants.all.2 plants.all.3 plants.all.4
  plants.neg plants.neg.0 plants.neg.1 plants.neg.2 plants.neg.3 plants.neg.4
  plants.neu plants.neu.0 plants.neu.1 plants.neu.2 plants.neu.3 plants.neu.4
  plants.pos plants.pos.0 plants.pos.1 plants.pos.2 plants.pos.3 plants.pos.4
  out constants
  fine_tbl
  large_tbl
]

turtles-own [
  genotype             ;; A list used to store alleles
  heterozygosity       ;; heterozygosity levels
  bins                 ;; for plotting purposes (to give colors according to heterozygosity levels)
  adaptation
  rd_coef
  maternal.gen.relatedness
  distance.travelled
 ]

plants-own [
 pollinators
 progeny_maternal
 progeny_paternal
]
seeds-own [            ;; seeds of both seed-types (no advantage to separating them as agents like with the mature plants)
  pollination.distance ;; To monitor how the pollination function works
  mother
  father
  siblings
]

links-own [
  link.relatedness
  link.length
  rd
]

patches-own [
  relatedness.patch
  original.value
  clear?
]

;; (profiler procedure removed for NetLogo Web — dev-only, requires the unsupported `profiler` extension)


  ;;#################################################################################
  ;;#################################################################################
  ;;######################                                     ######################
  ;;######################              Setup                  ######################
  ;;######################                                     ######################
  ;;#################################################################################
  ;;#################################################################################


to setup
  ca ;; clear the workspace
  reset-ticks
  ask patches [set pcolor 9.9]
  if (local-adaptations? or not (fragmentation = "no")) [landscape.setup]
  ask patches [hatch.seeds.setup ]

  set current-state "not-started"
  set-default-shape sides "line"
  set-default-shape seeds "circle"
  set-default-shape plants "flower"
  set-default-shape plants "flower"
  set selected no-turtles
  set patch_number count patches


end

to landscape.setup
    ask patches [
    if (random-float 1.0 < 0.5) [
        set pcolor 0]]
   repeat (10 * autocorrelation) [diffuse pcolor 0.25]


end

to hatch.seeds.setup
        sprout 1 [                                     ;; seeds that haven't randomly died produce one seed each:
          set breed seeds
    set size 0.25
          set genotype n-values 100 [shuffle [0 1] ]
          set heterozygosity random-float 1.0
          calc.bins
          set distance.travelled random-float 3.0
    set maternal.gen.relatedness random-float 1.0
          set bins 4
    if (local-adaptations?) [set adaptation [pcolor] of patch-here]
        ;set adaptation random-float 9.9
    ifelse (strategies = "single")
    [set rd_coef dispersal-relatedness]
    [set rd_coef (one-of read-from-string (word "[" multi-strategies "]"))]
  ]
end

to change.landscape
    if (ticks mod 5 = 0) [
  ask n-of (0.1 * patch_number) patches [set pcolor 9.9 ]
  ask n-of (0.1 * patch_number) patches [set pcolor 0]
  repeat 3 [diffuse pcolor autocorrelation]]
end



;; (commented-out desktop stats-table scaffolding removed for the NetLogo Web variant)


  ;;#################################################################################
  ;;#################################################################################
  ;;######################                                     ######################
  ;;######################          Main procedures            ######################
  ;;######################                                     ######################
  ;;#################################################################################
  ;;#################################################################################

to go
  kill.parents
  if (changing.landscape) [change.landscape]

  establish


  If (check.stop) [stop]

 if( ticks > 1 and plant_plots) [

       write_tables
    plots

       ]

 pollinate


 disperse.seeds

  tick

end


to establish ;; relatedness-depenedent establishment: more heterozygous seeds have higher establishment probabilities
   ask patches[
    if (any? seeds-here) [
      if (fragmentation != "no")  [ifelse (fragmentation = "black uninhabitable") [if (pcolor < 3) [ask seeds-here [die]]]
        [if (pcolor > 7) [ask seeds-here [die]]]]]

        if (any? seeds-here) [
           ifelse (local-adaptations?)
          [ifelse (coupled = true)
          [ask min-one-of seeds-here [abs ((heterozygosity * 10) - [pcolor] of patch-here)]
          [ifelse (abs ((heterozygosity * 10) - [pcolor] of patch-here) > 3) [die] [hatch.plants ]]]
          [ask min-one-of seeds-here [abs (adaptation - [pcolor] of patch-here)]
          [ifelse (abs (adaptation - [pcolor] of patch-here) > 3) [die] [hatch.plants ]]]]
      [ask one-of seeds-here [hatch.plants]]]]
  ask seeds [die]
  calc.fixation
end

to hatch.plants ;; procedure to create plants from seeds
  hatch 1 [
    set hidden? false
    set breed plants
    set size 1.5
    set progeny_maternal []
    set progeny_paternal []
    set pollinators []
  ]
end


;; Pure-NetLogo replacement for the weighted-n-of-with-repeats primitive (the `rnd` extension is unsupported on NetLogo Web).
;; pairs = list of [agent weight]; returns a list of k agents sampled WITH REPLACEMENT via roulette-wheel selection.
to-report wn-with-repeats [ pairs k ]
  let tot sum map [ p -> last p ] pairs
  let picks []
  if (empty? pairs or tot <= 0) [ report picks ] ;; guard empty / zero-weight cases so it never errors
  repeat k [
    let r random-float tot
    let cum 0
    let pick nobody
    foreach pairs [ p -> if pick = nobody [ set cum cum + last p  if r < cum [ set pick first p ] ] ]
    if pick = nobody [ set pick first last pairs ] ;; floating-point safety net: take the last candidate
    set picks lput pick picks
  ]
  report picks
end

to pollinate
  set mutations.overall []
  set pollination.distances []
  ask plants [
    let .neighbors other plants with [distance myself < neighborhood-size ];; make a list of potential donors in distance up to a value set in the main screen
    if (any? .neighbors) [
    ;; map runs in the focal (pollinating) plant's context, so `distance c` is the distance from the focal plant
    ;; to candidate donor c (distance is symmetric) — identical weighting to the desktop weighted sampler [ 1 / e ^ distance myself ]
    let pair-list map [ c -> (list c (1 / e ^ (distance c))) ] [self] of .neighbors ;; weight each candidate donor by 1 / e^(distance to the focal plant)
    foreach (wn-with-repeats pair-list (2 + random-poisson 6)) ;; choose number of seeds (random-poisson, mean + var = 6) and sample donors WITH REPLACEMENT by distance weight
                  [x -> hatch.seed.go x
                   set pollinators lput x pollinators]]
    ask (turtle-set progeny_maternal) [set siblings [progeny_maternal] of myself]]
  set pollination.distances [pollination.distance] of seeds
end

to hatch.seed.go [x]
  hatch-seeds 1 [
    set mother myself
    set father x
    ask mother [set progeny_maternal lput myself progeny_maternal]
    ask father [set progeny_paternal lput myself progeny_paternal]
    set hidden? true;; hide seeds
    set size 0.25
    if (local-adaptations?)
    [set adaptation one-of (list [adaptation] of mother  [adaptation] of father) ]
;    set adaptation random-normal (([adaptation] of .mother + [adaptation] of .father) / 2) mutations] ;; some mutation
    set genotype calc.genotype self father
    set maternal.gen.relatedness calc.relatedness self mother
    set heterozygosity calc.heterozygosity
    calc.bins
    set pollination.distance distance father ;;  *PLOTTING RELATED PROCEDURE*
    set rd_coef [rd_coef] of mother
    ;set rd_coef one-of (list [rd_coef] of mother [rd_coef] of father)
    set distance.travelled 0
    ]
end

to disperse.seeds

(ifelse
  (dispersal.strategy = "1")
    [ask plants [
        let .rd ifelse-value (strategies = "multi")
        [rd_coef]
        [dispersal-relatedness]
        (foreach sort-on [ .rd * maternal.gen.relatedness] (turtle-set progeny_maternal) sort n-values length progeny_maternal  [random-exponential (dist * 0.5)]
        [ [x y] -> ask x [
          set distance.travelled y]])]]
      (dispersal.strategy = "2") [
        let .min min [maternal.gen.relatedness] of seeds ;; to standardize dispersal-distances
        let .max max [maternal.gen.relatedness] of seeds ;; to standardize dispersal-distances
        let .range .max - .min ;; to standardize dispersal-distances
        let .sd standard-deviation [maternal.gen.relatedness] of seeds
    ask plants [
      ask (turtle-set progeny_maternal)
        [set distance.travelled (ifelse-value (strategies = "multi")
          [calc.dispersal.dist [rd_coef] of myself maternal.gen.relatedness .min .max .range .sd]
          [calc.dispersal.dist dispersal-relatedness maternal.gen.relatedness .min .max .range .sd])
    ]]])
      ask seeds [
        rt random 360 ;; choose a random dispersal direction
        forward distance.travelled]


end

to kill.parents
    ask plants [die]
 display
end


  ;;#################################################################################
  ;;#################################################################################
  ;;######################                                     ######################
  ;;######################       Reporters & Cakculators       ######################
  ;;######################                                     ######################
  ;;#################################################################################
  ;;#################################################################################

to-report calc.genotype [a b]

  let .genotype (map [ [ x y ] -> replace-item (random 2) x one-of y]
    [genotype] of a [genotype] of b)
  let .mutations round((random-exponential mutations )* 100)

  set mutations.overall lput .mutations mutations.overall
  repeat mutations [set .genotype replace-item random 100 .genotype one-of [ [0 1] [1 0] ]]
  report .genotype
end

to-report calc.heterozygosity
  report length (filter [x -> first x != last x] genotype) / 100
end

to calc.fixation
  set fixation.map []
  let denom (2 * count plants)
  set fixation.map map [x -> ifelse-value (all? plants [reduce - item x genotype = 0] and length remove-duplicates [item x genotype] of plants = 1)
    [0]
    [table:get-or-default (table:counts reduce sentence [item x genotype] of plants) 1 0 / denom]] (range 100)
end

to-report calc.relatedness [genotype1 genotype2]
  report (200 - sum (map [ [x y] -> (abs ( item 0 x + item 1 x - item 0 y - item 1 y))] [genotype] of genotype1 [genotype] of genotype2)) / 200
end


to-report calc.dispersal.dist [.B x .min .max .range .sd]
report (ifelse-value
    (.B = 1) [random-exponential ( dist * ((x - .min) / .range))]
      (.B = 0) [random-exponential (dist * (random-normal 0.5 .sd))]
        (.B = -1) [random-exponential ( dist * (( .max - maternal.gen.relatedness) / .range)) ] )
end

to calc.bins ;; use color codes to visualize dynamics
  set bins (ifelse-value
        (heterozygosity <= 0.2) [0] ;; most homozygous
        (heterozygosity <= 0.4) [1] ;;
        (heterozygosity <= 0.6) [2] ;;
        (heterozygosity <= 0.8) [3] ;;
         [4]) ;; most heterozygous
  set color (ifelse-value ;; color codes for figures
      (bins = 0) [red]
      (bins = 1) [orange]
      (bins = 2) [43]
      (bins = 3) [green]
       [blue])
end


to calc.sgs

  repeat 5 [
  let y random world-width
  let .sample plants-on patches with [pycor = y]
    ask .sample [
    create-links-with other .sample [
        set hidden? true
          set rd [rd_coef] of end1
        set link.relatedness calc.relatedness end1 end2]]]

end

to-report check.stop ;; Stop the simulation if population size reaches a critical thershold of 10%, if more than 90% of the loci are fixated, or if 2000 generations have passed
  let .stop (ifelse-value
  (ticks > 2000) [1]
  (count plants < (patch_number / 20)) [2]
  (strategies = "multi" and length remove-duplicates reduce sentence (list [rd_coef] of turtles) < 2) [3]
  [4] )
  let .output (ifelse-value
  (.stop = 1) ["Reached 2000 time steps"]
  (.stop = 2) ["Too few plants left"]
    (.stop = 3) [(word "Only " [rd_coef] of one-of plants " strategy left")]
  (.stop = 4) ["Stop conditions not fullfilled"])
  ifelse (.stop < 4) [set out lput .output out output-print .output report TRUE] [report FALSE]
end


  ;;#################################################################################
  ;;#################################################################################
  ;;######################                                     ######################
  ;;######################               Tables                ######################
  ;;######################                                     ######################
  ;;#################################################################################
  ;;#################################################################################


;; write_tables is a no-op in the NetLogo Web variant: it only ever populated the `stats`
;; extension tables that fed the regression/quantile analysis plots (omitted here).
;; The pure-NetLogo plots in `plots` do not depend on it.
to write_tables
end


  ;;#################################################################################
  ;;#################################################################################
  ;;######################                                     ######################
  ;;######################               Plots                 ######################
  ;;######################                                     ######################
  ;;#################################################################################
  ;;#################################################################################



to plots ;; plots
  let _plants []
  (ifelse   (plot? = "All") [
      set _plants plants]
[let y 0
      (ifelse
  (plot? = "Negative") [set y -1]
  (plot? = "Neutral") [set y 0]
  (plot? = "Positive") [set y 1])

      set _plants plants with [rd_coef = y]
    ask turtles [
        ifelse (rd_coef = y) [set hidden? false] [set hidden? true]]])


  ;;#########################################
  ;;#########################################

  ;; (fine-scale / large-scale SGS relatedness-by-distance regression plots and "sgs coefs dynamics"
  ;;  omitted in the web variant — they require the desktop `stats` extension)
  ;;#########################################
  ;;#########################################

  set-current-plot "Environmental variability"

  clear-plot
  set-plot-x-range 0 10
  set-histogram-num-bars 10
  set-current-plot-pen "default"
  histogram [pcolor] of patches

  ;;#########################################
  ;;#########################################

    set-current-plot "fixation"
    calc.fixation
    clear-plot
    set-current-plot-pen "pen-1" plot-pen-down plotxy 0 1 plotxy 100 1 plot-pen-up
    set-current-plot-pen "default" (foreach (range length fixation.map) fixation.map [ [x y] -> plotxy x y])
  ;;#########################################
  ;;#########################################

  set-current-plot ("fixation 2")
    plotxy ticks length filter [x -> x = 0 or x = 1] fixation.map



  ;;#########################################
  ;;#########################################


  set-current-plot ("mutation")
clear-plot
  set-plot-pen-interval 0.5
  set-histogram-num-bars 20
    histogram mutations.overall


  ;;#########################################
  ;;#########################################


  set-current-plot ("plants realized kernel" )

  clear-plot
  set-plot-y-range 0 count _plants
  set-plot-pen-interval 0.5
  set-histogram-num-bars 20
  histogram [distance.travelled] of _plants

  ;;#########################################
  ;;#########################################

   (foreach (list 4 3 2 1 0) (list blue green 43 orange red) [[x y] ->

  ;;#########################################
  ;;#########################################

      set-current-plot "plants quantities"

      set-plot-y-range 0 count _plants
      set-current-plot-pen  (word x)
        set-plot-pen-color y
          set-plot-pen-interval 0.5
  set-histogram-num-bars 20
    plotxy ticks count plants with [bins = x]

  ;;#########################################
  ;;#########################################

foreach [-1 0 1] [z ->

set-current-plot (word "plants mean dispersal distance coef = " z)
      let .plants plants with [rd_coef = z ]
let ..plants .plants with [bins = x]
if (any? ..plants) [
set-current-plot-pen (word x)
  set-plot-pen-color y
        plotxy ticks mean [distance.travelled] of ..plants

  set-current-plot "plants mean dispersal distance"

  set-current-plot-pen (word z)
              plotxy ticks  mean [distance.travelled] of .plants
      ]

  ;;#########################################
  ;;#########################################

      set-current-plot (word "plants realized kernel heterozygosity = " x)
      clear-plot
      set-plot-y-range 0 count .plants + 1

      set-plot-x-range 0 10
      set-histogram-num-bars 20
      set-plot-pen-interval 0.5
      set-plot-pen-color y
      histogram [distance.travelled] of ..plants
  ]])

  ;;#########################################
  ;;#########################################

  ;; ("Dispersal quantiles" plot omitted in the web variant — requires the `stats` extension)

  ;;#########################################
  ;;#########################################

  set-current-plot "plants heterozygosity"

  clear-plot
  set-plot-x-range 0 1
  set-plot-pen-interval 0.5
  set-histogram-num-bars 20
  histogram [heterozygosity] of _plants

  ;;#########################################
  ;;#########################################

  ;; ("heterozygosity dispersal plants" scatter + regression omitted in the web variant — requires the `stats` extension)

  ;;#########################################
  ;;#########################################

  ;; ("Heterozygosity median SE" quantile plot omitted in the web variant — requires the `stats` extension)

  ;;#########################################
  ;;#########################################

  ;; ("heterozygosity dispersal over time" regression-slope plot omitted in the web variant — requires the `stats` extension)

  ;;#########################################
  ;;#########################################



  ;; ("mother-offspring relatedness" scatter + regression omitted in the web variant — requires the `stats` extension;
  ;;  the "mother-offspring relatedness histogram" below is pure NetLogo and is kept)

  ;;#########################################
  ;;#########################################


set-current-plot "#plants by rd-strategy"

  (foreach [-1 0 1] ["-1" "0" "1"] [[x y] -> if (any? plants with [rd_coef = x]) [set-current-plot-pen y plotxy ticks count plants with [rd_coef = x]]])


  set-current-plot "mother-offspring relatedness histogram"
  clear-plot
  set-plot-x-range 0.5 1
  set-plot-pen-interval 0.5
  set-histogram-num-bars 20
  histogram [maternal.gen.relatedness] of _plants


  ;;#########################################
  ;;#########################################

   set-current-plot "pollination"
  histogram pollination.distances

end

to make-plot [.selected]
      set-current-plot "selected"

  (ifelse
  (what-to-plot = "sgs") [
     clear-plot
      let .list []
       set-current-plot-pen "x"
            ask .selected [
        create-links-with other .selected [
          set link.relatedness calc.relatedness end1 end2
          set link.length link-length
          plotxy link-length link.relatedness
          set .list lput (list link.relatedness link-length) .list
      ]]

    ;; (regression "mean" line for the sgs "selected" plot omitted — requires the `stats` extension;
    ;;  the pure-NetLogo relatedness-by-distance scatter above is kept)
      ask links [die]
    ]
  [set-plot-pen-mode 1
      set-plot-x-range 0 1
      set-plot-pen-interval 0.05
      (ifelse
  (what-to-plot = "mother-offspring") [histogram [maternal.gen.relatedness] of .selected]
      (what-to-plot = "adaptation") [histogram [adaptation / 10] of .selected]
      (what-to-plot = "heterozygosity") [histogram [heterozygosity] of .selected]
  (what-to-plot = "dispersal distances") [set-plot-x-range 0 10 set-plot-pen-interval 0.5 histogram [distance.travelled] of .selected])])

end

to plot_selected
  make-plot selected
end


  ;;#################################################################################
  ;;#################################################################################
  ;;######################                                     ######################
  ;;######################           Interactive               ######################
  ;;######################                                     ######################
  ;;#################################################################################
  ;;#################################################################################

to interactive_progeny
  if (mouse-down?) [
let .plants.here []
  let any.color random-float 140
  ask patch mouse-xcor mouse-ycor
  [ask plants-here
  [set color any.color
  set size 1.2
  create-links-with (turtle-set progeny_maternal) [
  set hidden? false
          set thickness 0.2
    set color [color] of end1]
    ask (turtle-set progeny_maternal)[
    set size 0.3
    set hidden? false
  ]
    set selected (turtle-set selected self progeny_maternal)]]]
display
end

to interactive_relatedness
  if (mouse-down?) [
let .here []
  let any.color random-float 140
  ask patch mouse-xcor mouse-ycor
  [ask one-of turtles-here
  [set size size * 1
  set .here other turtles in-radius 2
  let .levels map [x -> calc.relatedness x self * 100] sort .here
  ask .here [
  set size size * 1

  set color scale-color red (calc.relatedness self myself * 100) min .levels max .levels
]
  set selected (turtle-set selected .here self)]]]
display
end

to interactive_pollinators
  if (mouse-down?) [
let .plants.here []
  let any.color random-float 140
  ask patch mouse-xcor mouse-ycor
  [set .plants.here plants-here
  ask .plants.here [
  set color any.color
  set size 1.2

  create-links-with (turtle-set pollinators) [
  set hidden? false
    set color [color] of end1]
    ]]
    ask .plants.here [
    ask turtle-set pollinators[
    set color [color] of myself
    set size 1.3
    set hidden? false
]]
  set selected (turtle-set selected .plants.here [pollinators] of .plants.here)]
display
end

to rd_color
  ask turtles [(ifelse (rd_coef = 1) [set color blue]
(rd_coef = 0) [ set color grey]
(rd_coef = -1) [ set color red])]
end

to htz_color
  ask turtles [
     set color (ifelse-value ;; color codes for figures
      (bins = 0) [red]
      (bins = 1) [orange]
      (bins = 2) [43]
      (bins = 3) [green]
       [blue]) ]
end

to dispersal_landscape
  ask patches [set original.value pcolor]
let .global []
ask patches with [any? turtles-here][
set .global lput mean [distance.travelled] of
turtles-here .global]
set .global reduce sentence reduce sentence .global
ask patches [ifelse (any? turtles-here)
[set pcolor scale-color blue mean [distance.travelled] of
turtles-here max .global min .global ]
[set pcolor white]]
end

to relatedness_landscape
  ask patches [if (clear? != ticks) [
set original.value pcolor
set relatedness.patch "null"
ask patches with [any? turtles-here and any? turtles-on neighbors][
set relatedness.patch mean reduce sentence reduce sentence map
[x -> map [y -> calc.relatedness y x * 100] sort turtles-here]
sort (turtle-set turtles-here turtles-on neighbors)] ]]
ask patches [ ifelse (any? turtles-here)
[set pcolor scale-color blue relatedness.patch
max [relatedness.patch] of patches min [relatedness.patch] of patches ]
[set pcolor white]]

end

to reset_landscape
  ask patches [set pcolor original.value]
display
end



to interactive_fragment
  let .pro read-from-string user-input (word "Proportion of cells to turn " .color "?")
ask patches [if (random-float 1.0 < .pro) [
set pcolor read-from-string .color]]
end

to interactive_patches
  let populations read-from-string user-input "How many patches to make?"
let initial.pop.size read-from-string user-input "What size?"


  repeat populations [               ; number of blobs I want to  have
    let blob-maker nobody
    crt 1 [ set blob-maker self
       setxy random-xcor random-ycor]  ; set random position of "blob-makers"
    repeat initial.pop.size [               ; size of one blob (number of patches of the same color close one to another)
    ask blob-maker [
       ask min-one-of patches with-max [ pcolor ] [ distance myself ] [ set pcolor black]
  rt random 360
  fd 1
]
]
ask blob-maker [ die ]
]
end

to diffusing
  diffuse pcolor autocorrelation
end

to interactive_coloration
  if (mouse-down?) [
  ask patch mouse-xcor mouse-ycor
  [set pcolor read-from-string .color ]]
display
end

to hide [x]
  ask x[ set hidden? true]
end

to show_all
  set hidden? false
end

to clear_links
  clear-links
  deselect
       ask seeds [
       set size 0.25]
       ask plants [
       set size 0.5]

       set selected no-turtles
end

to clear_plot
  set-current-plot "selected" clear-plot

deselect
end




to go2
  ; we use a different procedure depending on which state we're in
  ; the procedures will update the `current-state` to the next value as the user clicks
  (ifelse
    (current-state = "not-started")     [ start-selecting ]
    (current-state = "selecting")       [ handle-selecting ]
  )
end

to start-selecting

  if mouse-down? [
    ; on the first click we see, we record the mouse position as the start of selection
    set select-x mouse-xcor
    set select-y mouse-ycor
    set current-state "selecting"
  ]
end

to handle-selecting
  ifelse mouse-down? [
    select select-x select-y mouse-xcor mouse-ycor]
  [ask sides [die]
    make-plot selected
    set current-state "not-started"]
  display
end

to deselect
  ask sides [ die ]
  ask selected [   set color (ifelse-value ;; color codes for figures
      (bins = 0) [red]
      (bins = 1) [orange]
      (bins = 2) [43]
      (bins = 3) [green]
       [blue]) ]
  set current-state "not-started"
end

to select [x1 y1 x2 y2]
  make-side x1 y1 x2 y1
  make-side x1 y1 x1 y2
  make-side x1 y2 x2 y2
  make-side x2 y1 x2 y2
  set selected (turtle-set selected turtles with [selected? xcor ycor] selected)
  ask selected [ set color red ]
end

to make-side [x1 y1 x2 y2]
  create-sides 1 [
    set color gray
    setxy (x1 + x2) / 2
          (y1 + y2) / 2
    facexy x1 y1
    set size 2 * distancexy x1 y1
  ]
end

to-report selected? [x y]
  let y-max max [ycor] of sides  ; largest ycor is where the top is
  let y-min min [ycor] of sides  ; smallest ycor is where the bottom is
  let x-max max [xcor] of sides  ; largest xcor is where the right side is
  let x-min min [xcor] of sides  ; smallest xcor is where the left side is
  ; report whether the input coordinates are within the rectangle
  report x >= x-min and x <= x-max and
         y >= y-min and y <= y-max
end


  ;;#################################################################################
  ;;#################################################################################
  ;;######################                                     ######################
  ;;######################           Output data               ######################
  ;;######################                                     ######################
  ;;#################################################################################
  ;;#################################################################################


;; (dead reporters `coefs` and `quantiles` removed — they depended on the `stats` extension
;;  and were only referenced by the omitted analysis plots / BehaviorSpace metrics)





;  let .count.total count .plants
;  foreach [0 1 2 3 4] (list .bins.0 .bins.1 .bins.2 .bins.3 .bins.4) [ y z -> if (any? .plants with [bins = y) [let z .plants with [bins
;  let progeny_size count .progeny
;  let progeny_mean_dispersal mean [distance.travelled] of .progeny
;  let progeny_mean_M-O_relatedness mean [maternal.gen.relatedness] of .progeny
;  let progeny_mean_heterozygosity mean [heterozygosity] of .progeny
;  let radius1 mean map [x -> calc.relatedness self x] ([self] of plants in-radius 1 )
;  let radius2 mean map [x -> calc.relatedness self x] ([self] of plants in-radius 2 )
;  let radius3 mean map [x -> calc.relatedness self x] ([self] of plants in-radius 3 )
;  let radius4 mean map [x -> calc.relatedness self x] ([self] of plants in-radius 4 )
;  let radius5 mean map [x -> calc.relatedness self x] ([self] of plants in-radius 5 )
;  let fathers.relatedness mean map [x -> calc.relatedness self x] .fathers
;  let .tmp []
;  ask .progeny [ask other .progeny [set .tmp lput calc.relatedness self myself .tmp]]
;  let among-progeny-relatedness mean .tmp
;  let patch_color pcolor
;  let pollinators_mean_distance mean map [x -> distance x] .fathers
;  let count_pollinators sum sort remove-duplicates .fathers
;    [omitted-analytics add] plants_data (list distance.travelled progeny_mean_heterozygosity maternal.gen.relatedness xcor ycor
;    progeny_size progeny_mean_dispersal among-progeny-relatedness progeny_mean_heterozygosity radius1
;    radius2 radius3 radius4 radius5 rd_coef adaptation patch_color fathers.relatedness progeny_mean_M-O_relatedness count_pollinators pollinators_mean_distance)
;end

;to output.setup
;  let progeny_mean_heterozygosity []
;  let progeny_size []
;  let progeny_mean_dispersal []
;  let among-progeny-relatedness []
;  let radius1 []
;  let radius2 []
;  let radius3 []
;  let radius4 []
;  let radius5 []
;  let fathers.relatedness []
;  let progeny_mean_M-O_relatedness []
;  let count_pollinators []
;  let pollinators_mean_distance []
;end
;
;to output.data [progeny .fathers]
;  set progeny_size lput count progeny progeny_size
;  set progeny_mean_dispersal lput mean [distance.travelled] of progeny progeny_mean_dispersal
;  set progeny_mean_M-O_relatedness lput mean [maternal.gen.relatedness] of progeny progeny_mean_M-O_relatedness
;  set progeny_mean_heterozygosity lput mean [heterozygosity] of progeny progeny_mean_heterozygosity
;  set radius1 lput mean map [x -> calc.relatedness self x] ([self] of plants in-radius 1 ) radius1
;  set radius2 lput mean map [x -> calc.relatedness self x] ([self] of plants in-radius 2 )  radius2
;  set radius3 lput mean map [x -> calc.relatedness self x] ([self] of plants in-radius 3 ) radius3
;  set radius4 lput mean map [x -> calc.relatedness self x] ([self] of plants in-radius 4 ) radius4
;  set radius5 lput mean map [x -> calc.relatedness self x] ([self] of plants in-radius 5 ) radius5
;  set fathers.calc.relatedness lput mean map [x -> calc.relatedness self x] .fathers fathers.relatedness
;  set .tmp []
;  ask progeny [ask other progeny [set .tmp lput calc.relatedness self myself .tmp]]
;  set among-progeny-relatedness lput mean .tmp among-progeny-relatedness
;  set patch_color lput pcolor patch_color
;  set pollinators_mean_distance lput mean map [x -> distance x] .fathers pollinators_mean_distance
;  set count_pollinators lput sum sort remove-duplicates .fathers count_pollinators
;end
;
;to-report output.means .list
;  report mean (reduce sentence .list)
;end
@#$#@#$#@
TEXTBOX
255
3
995
33
Spatial genetic structure under habitat fragmentation. Plant colour = genetic diversity (heterozygosity):  red = inbred / homozygous  ->  blue = outbred / heterozygous. As fragmentation limits seed dispersal, related & inbred plants cluster in space.
11
0.0
1

GRAPHICS-WINDOW
255
35
415
197
-1
-1
13.71
1
10
1
1
1
0
1
1
1
-5
5
-5
5
1
1
1
ticks
30.0

BUTTON
0
400
56
435
setup
setup\n
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
58
400
114
434
NIL
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

BUTTON
114
400
170
435
go once
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
59
436
115
470
pollinate
pollinate
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
117
433
172
470
disperse seeds
disperse.seeds
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
1374
133
1688
253
plants mean dispersal distance
ticks
dispersal distance
0.0
1.0
0.0
1.0
true
true
"" ""
PENS
"0" 1.0 0 -7500403 true "" ""
"-1" 1.0 0 -2674135 true "" ""
"1" 1.0 0 -13345367 true "" ""

PLOT
894
293
1054
413
plants heterozygosity
heterozygosity
# plants
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 1 -8630108 true "" ""

SWITCH
890
10
980
43
plant_plots
plant_plots
0
1
-1000

BUTTON
1070
10
1126
44
plot
plots\n
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
0
433
56
468
establish
  kill.parents \n  if (changing.landscape) [change.landscape]\n  establish\n   if (ticks > 1)[ write_tables ]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
1054
173
1214
293
pollination
Distance
Quantity
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 1 -16777216 true "" ""

PLOT
733
173
893
293
plants quantities
ticks
# plants
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"0" 1.0 0 -7500403 true "" ""
"1" 1.0 0 -2674135 true "" ""
"2" 1.0 0 -955883 true "" ""
"3" 1.0 0 -6459832 true "" ""
"4" 1.0 0 -1184463 true "" ""

BUTTON
170
398
226
432
NIL
clear-all
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SWITCH
190
563
280
596
local-adaptations?
local-adaptations?
1
1
-1000

PLOT
1054
293
1214
413
Environmental variability
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 1 -7171555 true "" ""
"rd" 1.0 1 -8630108 true "" ""
"ri" 1.0 1 -5825686 true "" ""

PLOT
1214
13
1374
133
plants realized kernel
distance travelled
# seeds
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 1 -8630108 true "" ""
"mean" 1.0 1 -16777216 true "" ""

PLOT
1214
133
1374
253
plants realized kernel heterozygosity = 0
distance travelled
# seeds
0.0
20.0
0.0
10.0
true
false
"" ""
PENS
"mean" 1.0 1 -2674135 true "" ""

PLOT
1214
253
1374
373
plants realized kernel heterozygosity = 1
distance travelled
# seeds
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"mean" 1.0 1 -955883 true "" ""

PLOT
1214
373
1374
493
plants realized kernel heterozygosity = 2
distance travelled
# seeds
0.0
20.0
0.0
10.0
true
false
"" ""
PENS
"mean" 1.0 1 -4079321 true "" ""

PLOT
1214
493
1374
613
plants realized kernel heterozygosity = 3
distance travelled
# seeds
0.0
20.0
0.0
10.0
true
false
"" ""
PENS
"mean" 1.0 1 -10899396 true "" ""

PLOT
1214
613
1374
733
plants realized kernel heterozygosity = 4
distance travelled
# seeds
0.0
20.0
0.0
10.0
true
false
"" ""
PENS
"mean" 1.0 1 -13345367 true "" ""

TEXTBOX
1278
38
1333
59
Total
20
0.0
1

TEXTBOX
1273
158
1349
180
HtZ = 0
20
15.0
1

TEXTBOX
1268
519
1345
541
HtZ = 3
20
55.0
1

TEXTBOX
1268
394
1353
416
HtZ = 2
20
44.0
1

TEXTBOX
1268
634
1351
656
HtZ = 4
20
105.0
1

TEXTBOX
1268
279
1342
301
HtZ = 1
20
25.0
1

SLIDER
95
503
188
536
neighborhood-size
neighborhood-size
1
10
3.0
1
1
NIL
HORIZONTAL

PLOT
1053
533
1213
653
fixation
NIL
NIL
0.0
100.0
0.0
1.1
true
false
"" ""
PENS
"default" 1.0 1 -16777216 true "" ""
"pen-1" 1.0 0 -2674135 true "" ""

PLOT
1054
413
1214
533
fixation 2
NIL
NIL
0.0
10.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" ""

SLIDER
285
518
378
551
autocorrelation
autocorrelation
0
1
0.2
0.1
1
NIL
HORIZONTAL

SLIDER
94
471
187
504
mutations
mutations
0
1
0.0
0.01
1
NIL
HORIZONTAL

PLOT
733
293
893
413
mutation
NIL
NIL
0.0
99.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 1 -16777216 true "" ""

SWITCH
290
553
380
586
changing.landscape
changing.landscape
1
1
-1000

SWITCH
190
598
280
631
Coupled
Coupled
1
1
-1000

CHOOSER
190
518
283
563
fragmentation
fragmentation
"no" "black uninhabitable" "white uninhabitable"
0

PLOT
893
173
1053
293
mother-offspring relatedness histogram
NIL
NIL
0.5
1.0
0.0
10.0
true
false
"" ""
PENS
"default" 0.05 1 -16777216 true "" ""

OUTPUT
1742
115
2624
603
12

PLOT
1374
13
1688
133
#plants by rd-strategy
time
# plants
-2.0
2.0
0.0
10.0
true
true
"" ""
PENS
"-1" 1.0 0 -2674135 true "" ""
"0" 1.0 0 -7500403 true "" ""
"1" 1.0 0 -13345367 true "" ""

SLIDER
800
10
893
43
sgs.scale
sgs.scale
1
9
3.5
0.5
1
NIL
HORIZONTAL

BUTTON
735
10
800
44
reset y axis
set-current-plot \"plants mean dispersal distance\"\nset-plot-y-range 0 1\n\nset-current-plot \"plants mean dispersal distance coef = -1\"\nset-plot-y-range 0 1\n\nset-current-plot \"plants mean dispersal distance coef = 0\"\nset-plot-y-range 0 1\n\nset-current-plot \"plants mean dispersal distance coef = 1\"\nset-plot-y-range 0 1\n
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
0
653
92
686
switch_p
switch_p
0
1
0.1
0.1
1
NIL
HORIZONTAL

CHOOSER
0
575
93
620
strategy
strategy
"positive" "neutral" "negative"
2

BUTTON
95
653
151
687
switch
ask n-of (switch_p * count turtles) turtles \n[set rd_coef (ifelse-value \n(strategy = \"positive\") [1]\n(strategy = \"neutral\") [0]\n(strategy = \"negative\") [-1]\n)]
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
95
538
188
571
dist
dist
0
10
4.0
0.5
1
NIL
HORIZONTAL

SLIDER
0
618
93
651
dispersal-relatedness
dispersal-relatedness
-1
1
1.0
1
1
NIL
HORIZONTAL

PLOT
1374
253
1689
373
plants mean dispersal distance coef = 1
time
Dispersal distance
0.0
10.0
0.0
1.0
true
false
"" ""
PENS
"0" 1.0 0 -16777216 true "" ""
"1" 1.0 0 -7500403 true "" ""
"2" 1.0 0 -2674135 true "" ""
"3" 1.0 0 -955883 true "" ""
"4" 1.0 0 -6459832 true "" ""

PLOT
1374
373
1691
493
plants mean dispersal distance coef = 0
NIL
NIL
0.0
10.0
0.0
1.0
true
false
"" ""
PENS
"0" 1.0 0 -16777216 true "" ""
"1" 1.0 0 -7500403 true "" ""
"2" 1.0 0 -2674135 true "" ""
"3" 1.0 0 -955883 true "" ""
"4" 1.0 0 -6459832 true "" ""

PLOT
1374
493
1691
613
plants mean dispersal distance coef = -1
NIL
NIL
0.0
10.0
0.0
1.0
true
false
"" ""
PENS
"0" 1.0 0 -16777216 true "" ""
"1" 1.0 0 -7500403 true "" ""
"2" 1.0 0 -2674135 true "" ""
"3" 1.0 0 -955883 true "" ""
"4" 1.0 0 -6459832 true "" ""

CHOOSER
0
531
93
576
strategies
strategies
"single" "multi"
0

INPUTBOX
0
471
93
531
multi-strategies
-1 0 1
1
0
String

BUTTON
170
258
225
292
select
go2
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
0
13
234
177
selected
NIL
NIL
0.0
1.0
0.5
1.0
true
false
"" ""
PENS
"x" 1.0 2 -8630108 true "" ""
"mean" 1.0 0 -16777216 true "" ""

BUTTON
58
226
114
260
clear plot
clear_plot
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

CHOOSER
0
180
93
225
what-to-plot
what-to-plot
"sgs" "mother-offspring" "adaptation" "heterozygosity" "dispersal distances" "relatedness"
0

BUTTON
0
226
56
259
plot selected
plot_selected
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
58
260
114
294
progeny
interactive_progeny\n
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
114
226
197
259
clear links
clear_links
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
113
258
169
292
Pollinators
interactive_pollinators
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
0
327
56
361
rd color
rd_color
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
58
327
114
361
htz color
htz_color
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
0
362
56
396
disperal
dispersal_landscape
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
0
260
56
294
relatedness
interactive_relatedness
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
114
295
170
329
show all
show_all
NIL
1
T
TURTLE
NIL
NIL
NIL
NIL
1

BUTTON
58
295
114
329
hide seeds
hide seeds
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
0
293
56
327
hide plants
hide plants
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
58
363
114
397
relatedness
relatedness_landscape
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
114
328
170
362
diffuse
diffusing
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
380
545
436
579
fragment
interactive_fragment
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
380
510
437
544
patches
interactive_patches
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
380
475
435
510
manual coloration
interactive_coloration
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

CHOOSER
285
473
381
518
.color
.color
"black" "white"
1

BUTTON
113
362
169
396
reset
reset_landscape
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
175
330
230
390
plot.every
1.0
1
0
Number

CHOOSER
190
473
283
518
dispersal.strategy
dispersal.strategy
"1" "2"
0

CHOOSER
980
10
1072
55
plot?
plot?
"All" "Negative" "Neutral" "Positive"
0

@#$#@#$#@
# ⚠️ Browser (NetLogo Web) demo — REDUCED version

**This is a reduced, browser-compatible LITE demo of the desktop model.** The simulation dynamics are **identical** to the full desktop version (`habitat-fragmentation-dispersal-model.nlogo`, in the same repository): the same genetics, pollination, dispersal, fixation and mutation engine runs here unchanged. The weighted pollen-donor sampling that the desktop model performs with the `rnd` extension is reproduced here with an equivalent pure-NetLogo roulette-wheel helper (`wn-with-repeats`), so results are statistically the same.

What is **omitted** in this browser demo (because NetLogo Web does not support the `stats`, `rnd` and `profiler` extensions): the 8 regression / quantile **analysis plots** that relied on the desktop `stats` extension — *fine-scale* and *large-scale* SGS regressions, *sgs coefs dynamics*, *Dispersal quantiles*, *Heterozygosity median SE*, *heterozygosity dispersal plants*, *heterozygosity dispersal over time*, and the *mother-offspring relatedness* scatter+regression. All other plots (realized dispersal kernels, fixation, mutation, heterozygosity histograms, plant quantities, mean dispersal distances, pollination, environmental variability, etc.) run normally. BehaviorSpace experiments and the developer profiling button are also removed (the browser has no BehaviorSpace).

➡️ **For the complete model with the full statistical-analysis plots and BehaviorSpace experiments, open the desktop version `habitat-fragmentation-dispersal-model.nlogo` in NetLogo desktop.**

> Tip: the `plant_plots` switch is ON by default so the surviving plots draw as the model runs.

---

# My Model Description


## 1. Purpose and patterns
_The aim_ of of the model is to study the long-term effects and dynamics of the effects of relatedness on dispersal-strategies, and the potential feedback loop with the population’s spatial genetic structure. Specifically, the aim of this version of the model is to see if, and under which conditions, a relatedness-dependent dispersal strategy is more advantageous than a relatedness-INdependent dispersal strategy.

**Background:** The model is ment to elaborate on results of a common-garden experiment, in which we have conducted manual cross-pollinations using pollen-donors of increasing geographic distances (self, kin, near neighbors, distant neighbors, etc.), working under the assumption that geographic distance is a proxy of genetic distance (isolation by distance). In this work we found that as the geographic (and probably genetic) distance between the parent-plants increases, there is an increase in the mass of the produced seeds, and thus a decrease in their dispersal potential (because these are wind-dispersed seeds, so the heavier they are the closer-by they fall). Theoretically, this may be caused by an inbreeding load that causes seed mass reduction.

**Theoretical background:** Theoretically, a relatedness-dependent dispersal strategy could be advantagous. I am basing this hypothesis on a few assumptions:
1. In plants, often both seed- and pollen-dispersal are distance-dependent.
2. Dispersal is maternally-controlled, as the dispersal-related structures of the seed are constructed froma a maternal tissue sorrounding the embryo.
3. multi-paternity is ubiquitous in plants - often there are multiple pollen-donors within the same progeny/capitulum. Thus, we may expect heterogeneous levels of genetic-relatedness to exist between the mother-plant and the seeds, as well as among the seeds.

In general, theory predicts dispersal to evolve as a balance between the costs and benefits of dispersal, which are predicted to vary depending on which selection pressures are acting in the system. I suggest that these costs and benefits may be relatedness-dependent, and thus, theoretically, it may be advantageous to have a flexible, relatedness-dependent, dispersal strategy. For example: if kin competition is a strong selection pressure, it may be more advantageous to have the seeds that are more genetically-similar disperse further away, to reduce levels of competition among them and to maximize inclusive-fitness. On the other hand, if local adaptations are a strong selection pressure, we may expect an opposite pattern, as it will be more advantageous to have the seeds that are more genetically-similar to the mother-plant stay closer by, as they are probably more well-adapted to the environmental conditions in the maternal habitat.


#### Research questions: 
  * Is a relatedness-dependent dispersal strategy superior to a relatedness-independent one?
  * Does a relatedness-dependent dispersal strategy contribute to the maintenence of low inbreeding levels?
  * What is the spatial-genetic struture of a population that has a relatedness-dependent dispersal strategy?

#### Predictions:

  * If relatedness between parent-plants affects dispersal distances, and since pollination is distance dependent, overall inbreeding levels and levels of kin competition will be reduced, and the spatial genetic structure will be weak.
    * _To test this prediction:_ Run simulations with and without effect of relatedness on dispersal distance
    * Compare spatial genetic structures after X generations
    * compare fitness of relatedness-dependent and relatedness-independent dispersal strategies (proportion of established offspring)

  * A relatedness-dependent dispersal strategy will be superior to a relatedness-independent dispersal strategy, due to the negative effects of inbreeding depression
    * _To test this prediction:_ Run simulations where individuals that have a relatedness-dependent dispersal strategy “compete” with individuals who have a relatedness-independent dispersal strategy 
    * see which strategy takes over 
    * Maybe not one will take over, but both will co-exist. In such case, see what long-term dynamics emerge between the two strategies. Do we reach a stable state?

## 2. Entities, State Variables, and Scales
 

There are 5 entities in the model: patches, plants, pollen-grains, female gametes and seeds.
Plants are annual, wind-dispersed and wind-pollinated, outcrossed-only.

#### State variables:
•	Patches – location
•	Seeds - home patch (the patch of the mother plant), current patch (the patch to which the seed was dispersed), mass, diploid genome (9 diallelic genes)
•	Plants – have the same variables as seeds
•	Pollen grains - parent plant, haploid genome and an indicator for whether or not it was already used for pollination
•	Female gametes - parent plant and a haploid genome.
In addition, all entities have an indicator stating whethe their dispersal-strategy is relatedness-dependent (rd) or relatedness-independent (ri). 

#### Scales 
One time-step equals to one life-time unit. I am not sure about the extent at this point. 
One grid cell represents a habitat suitable for 1 individual. Again, unsure about the extent.




## 3. Process Overview and Scheduling

In each time step: 

  1. Establishment: In each cell one seed establishes. Establishment probability is proportional to the masses of the seeds. After establishment all seeds die, and a certain percentage of established plants are randomly chosen and die as well (this proportion can be manipulated using a slider).

  2. Production of reproduction units: Plants produce female gametes and pollen grains. The number of gametes and grains produced by each individual is drawn from a Poisson-distribution with mean that can be manipulated using a slider. Each gamete and grain inherit half of their parent’s genotype, so they are haploid. 

  3. Pollination: Pollination can only occur between plants that have the same dispersal-strategy. Each plant collects pollen from other plants in a distance-dependent manner. The number of potential donors and the neighborhood-size (max distance of allowable donors) can be manipulated using sliders). Then, to produce a seed, each gamete of the plant randomly chooses one of the collected grains, after checking it has not been previously used for pollination. (Actually, I need to recode this part, because I think currently one pollen grain can be used by several plants).

  4. Seed production: The seed inherits one set of alleles from the gamete and one from the chosen pollen-grain, thus it has a diploid genome. Then, an inbreeding coefficient is calculated as a function of homozygosity levels (non-linear relationship between homozygosity and inbreeding. The exact shape of the function can be manipulated with a slider). Then, seed mass is given by drawing a random number from a normal distribution. The mean value of this distribution depends on the dispersal strategy: 
For the relatedness-dependent strategy, the mean value is the value of homozygosity
For the relatedness-independent strategy, the mean value is drawn from a poisson distribution whose mean value is maximum-homozygosity/2 
(so that in both cases the mean mass will be the same)
Finally, for plotting purposes, the plants are given colors as a function of inbreeding levels.

  4. Dispersal: Each seed disperses in a random direction and to a distance that is a function of its’ mass.

  5. After seed dispersal, all mature plants in the landscape die

  6. Spatial genetic structure estimation: each plant calculates its’ genetic distance from all other plants who share the same dispersal-strategy by comparing their genomes and seeing in how many microsatellites they are identical. Then geographic distance is measured, and a plot is made (but for only one individual at this point, I think).


## 4. Design Concepts
 

## 5. Initialization
 
X Seeds are created and randomly dispersed in the landscape. Individuals are assigned random genotypes. Seeds turn to plants (or die) in cells (one per cell maximum), according to their mass relative to the mass of all other seeds within the cell (heaviest seed has the higest probability to establish).

## 6. Input Data
 
## 7. Submodels

#### population setup

Create N seeds
Place each seed in a random patch
Set each seed’s home-patch as the current patch it’s on 
Assign half of the seeds to a relatedness-dependent strategy, and the rest to a relatedness-independent strategy 
Assign each microsattelite with a value of 0, 1 or 2 ;; The value of each microsatelite (ms) represents the sum of values of the two alleles within it, where each allele can have a value of either 0 or 1.
Set homozygosity & mass ;; explained in the calc-homozygosity procedure

#### Calc homozygosity
;; To set homozygosity of seed based on its' genotype:
Set base homozygosity value as 1
Check the value of each microsatellite. If it is 0 or 2 -> add 1 to homozygosity
Repeat for all microsatellites.

;; To set seed mass:
Draw from a normal distribution with a mean of: 
if rd = true -> X * e ^ ( - random-poisson (homozygosity) / q )
if rd = false -> X * e ^ ( - random-poisson (5.5) / q )
;; where X is a value for seed mass that can be changed in the main screen, and q is a constant altering the decline steepness of the relationship-function, which can be changed in the main screen 

;; To set inbreeding-load:
 If 1 <= homozygosity <= 2 set inbreeding-load 0
 If 3 < homozygosity <= 4 set inbreeding-load 1
 If 5 < homozygosity <= 6 set inbreeding-load 2
 If 7 < homozygosity <= 8 set inbreeding-load 3
 If 9 < homozygosity <= 10 set inbreeding-load 4


#### establish

Each patch that has seeds on it chooses one seed. The probability of a seed being chosen is proportional to its mass.
The chosen seed hatches one plant, and dies.
The rest of the seeds on the patch die.
A proportion of X plants in the landscape are randomly chosen and killed. 
Patches with a plant on them are colored green. Patches without a plant on them are colored brown.


#### Produce gametes

Each plant produces X gametes. X is a number randomly drawn from a Poisson distribution with a mean of gametes-number. ;; gametes-number can be set in the main screen
give genotypes to gametes ;; explained in the give-genotype procedure
Set parent as the mother-plant


#### Produce pollen

Each plant produces X pollen grains. X is a number randomly drawn from a Poisson distribution with a mean of gametes-number * 2.
give genotypes to pollens
Set parent as the mother-plant
Set pollinated? False ;; used to account for whether the grain was already used for pollination or not

#### Give genotype
Each gamete/pollen inherits half of its parents alleles in the following manner:
If ms-a of the mother-plant = 2, set allele-a of gamete to 1 
If ms-a of the mother-plant = 1,  set allele-a of gamete to either 1 or 0 randomly 
Otherwise, set allele-a of gamete to 0
Repeat for allele-b to allele-i


#### Produce seeds
;; To make a list of candidate pollen-grains, chosen in a distance-dependent manner:
Ask plants to make a list of patches made of:
 X immediate neighbors that have on them other plants with the same rd value as the asking plant
 X / 2 of patches found in a radius larger than 2 and smaller than Y / 2 that have on them other plants with the same rd value as the asking plant
 X / 4 of patches found in a radius larger than Y / 2 and smaller than Y that have on them other plants with the same rd value as the asking plant
;; X is the number of donors and Y is neighborhood size. Both can be set in the main screen
Ask plants to make a list of pollens from pollen-grains found on patches that are in the patch list, and who have the same rd value as the asking plant
Each plant asks each of its gametes to take a random pollen-grain from the list of pollen grains, with pollinated? = false & parent-of-pollen != parent-of-gamete, and use it to create a seed.
If there is no pollen grain that fulfills these conditions, the gamete dies
To create the seed’s genotype, sum the value of alleles of the pollen and the gamete in each locus
Set the pollen’s pollinated? to True


#### To disperse seeds
Ask pollens to die
Ask seeds to set home-patch to the current patch they’re on
If mass > 0 
 turn to a random direction
 Move forward X / mass cells. X is a number drawn from a poisson distribution with a mean of dispersal-distance
Set distance-travelled as the distance from current patch to home-patch


#### To kill parents

ask plants to die


#### To calculate spatial genetic structure
 
Ask each plant to make a list of all other plants with the same rd? value as itself
For each plant on the list, calculate genetic distance:
 If |ms-a-of-other-plant – ms-a-of-myself| = 2, set genetic-distance + 1/9
 If |ms-a-of-other-plant – ms-a-of-myself| = 1, set genetic-distance + 0.5/9
 Repeat for ms-b to ms-i
Calculate geographic-distance as the eucalidean distance between the two plants
Put both values in a list
Repeat for all plants
 

## 8. Things I need help with

  * Is there a more proper way to estimate relatedness? Like creating a pedigree of some sort?
    * Currently I am modelling relatedness by looking at homozygosity levels, but this is not really accurate, at least not with such few alleles. I could theoretically make 1,000 alleles instead of "just" 9, but I guess that would be computationally heavy.
    * Also, since I am using homozigosity levels for two purposes (calculating relatedness and assigning inbreeding-load which causes mass reduction), and this mix might have some artifacts I am not currently aware of.
  * Need to improve setup of multiple populations
    * The current way in which I am setting up the populations in the initialization phase is causing errors, especially when the landscape is small
  * How to output spatial genetic structure of all individuals, and not just one plant?
    * Currently I am only able to produce a plot for 1 plant. What I really want is to have a plot of mean sgs for all plants.
  * Improving speed of run?
    * Seems that right now I can only work on very small scales. Maybe I can increase the models' efficiency so I can work on larger landscapes?
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
NetLogo 6.2.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
</experiments>
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
1
@#$#@#$#@
