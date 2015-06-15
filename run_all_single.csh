#!/bin/csh

set BRlist="0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0"

foreach BR(${BRlist})

./prepare_draw_macro_for_single.csh $BR

end
