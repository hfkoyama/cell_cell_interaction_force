#Shellscript for data assimilation
#Written by Hiroshi Koyama (National Institute for Basic Biology, Japan)


#User should set time frames to be analyzed.
start_t=0
end_t=0

#file name for run
run_name='./run_data_assimilation_221211_06'

#optional: analysis ID
id=1

#run though time frame from start_t to end_t
for t in `seq $start_t $end_t`
do
	echo time_frame=${t}
#	run="${run_name} ${t}"
	run="${run_name} ${t} ${id}"
	${run}
done


