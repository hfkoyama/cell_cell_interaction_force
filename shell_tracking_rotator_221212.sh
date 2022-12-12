#Shellscript for tracking data rotator
#Written by Hiroshi Koyama (National Institute for Basic Biology, Japan)


#User should set time frames to be analyzed.
start_t=1
end_t=2

#file name for run
run_name='./run_tracking_rotator_221212_02'

#run though time frame from start_t to end_t
for t in `seq $start_t $end_t`
do
	echo time_frame=${t}
	run="${run_name} ${t}"
	${run}
done


