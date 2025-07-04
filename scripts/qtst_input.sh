#!/bin/bash
#---------------------------------------
# Generate the qtst analysis input file.
#---------------------------------------
cat > qtst_input.yaml <<!
deltaE0: 0.0486
i_omega: 17.092786
#down_limit: -100
temperatures: [200, 300, 400, 500, 600, 800, 1000, 1400, 1800, 2200]
#output_file: qtst_tunneling.yaml
!