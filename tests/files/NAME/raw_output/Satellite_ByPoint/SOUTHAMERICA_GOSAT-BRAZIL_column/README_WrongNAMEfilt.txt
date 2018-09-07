Meant to only level_order NAME filter to input data to remove entries where second GOSAT level had
a higher pressure than the NAME surface pressure.
Actually was applying all three filters (cutoff,level_order,dpressure_range) with
cutoff <= 5.0,layer_range betweem 50 - 500 m.
