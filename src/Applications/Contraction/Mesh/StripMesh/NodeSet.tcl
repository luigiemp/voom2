# Get Nodal position for selected nodes

*createmarkpanel nodes 1 "Select Nodes"
set myNodes [ hm_getmark nodes 1 ]
set fp [ open "NodeSet.txt" "w" ]

set nNodes [ llength $myNodes ]
puts "Number of Nodes: $nNodes"

foreach id $myNodes {
	puts $fp $id
}

close $fp
*clearmark nodes 1
