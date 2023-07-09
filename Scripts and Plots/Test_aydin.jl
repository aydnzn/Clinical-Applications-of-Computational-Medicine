############ open and close actibelt
using ActIO.Usb
act = open(find_rct2())
state(act)
datamode!(act)
recording_practice!(act, StartNoStop)
recording_practice!(act, StopAny)
state(act)

using Plots
using ActIO.File
using ActAnalysis.StepSlc
using DataFrames
using CSV
using ActAnalysis.Activity
using ActAnalysis.SpeedSvr
using ActAnalysis.Stepwave
using ActAnalysis.StepSlc


r = recording("oussama.hdr")
acc = r[accel,30000:120000]
acc_1 = acc[:,1]
acc_2 = acc[:,2]
acc_3 = acc[:,3]


stepslcdetector = StepSlcDetector()
slcsteps = stepslcdetector(acc)
DataFrame(slcsteps)
stepanalyser = StepAnalyser()
steps = stepanalyser(acc)
DataFrame(steps)

group_steps = steps[985:end,:]
group_steps_right = group_steps[1:2:end,:]
group_steps_left = group_steps[2:2:end,:]
size_group = 1:402;
#counter= 0;
MAT_acc_step_1=zeros(200,1)
MAT_acc_step_2=zeros(200,1)
MAT_acc_step_3=zeros(200,1)
MAT_acc_step_1_left=zeros(200,1)
MAT_acc_step_2_left=zeros(200,1)
MAT_acc_step_3_left=zeros(200,1)

for i in size_group

        start_index_of_step = datetime_to_index(group_steps.time[i],r.start,10) ;
        end_index_of_step = datetime_to_index(group_steps.time[i+1],r.start,10) -1;
        
        acc_step = r[accel,start_index_of_step:end_index_of_step];
        acc_step_1 = acc_step[:,1];
        acc_step_2 = acc_step[:,2];
        acc_step_3 = acc_step[:,3];
        
        size_acc_step_1 =size(acc_step_1,1);
        size_acc_step_2 =size(acc_step_2,1) ;
        size_acc_step_3 =size(acc_step_3,1);
        
        
        itp = LinearInterpolation(1:size_acc_step_1, acc_step_1) ;# create interpolation function
        range_acc_step_1 = range(1,stop=size_acc_step_1,length=200);
        interpolated_acc_step_1 = itp(range_acc_step_1);
        MAT_acc_step_1 = hcat(MAT_acc_step_1,interpolated_acc_step_1);
        
        itp = LinearInterpolation(1:size_acc_step_2, acc_step_2); # create interpolation function
        range_acc_step_2 = range(1,stop=size_acc_step_2,length=200);
        interpolated_acc_step_2 = itp(range_acc_step_2);
        MAT_acc_step_2 = hcat(MAT_acc_step_2,interpolated_acc_step_2);
        
        itp = LinearInterpolation(1:size_acc_step_3, acc_step_3) ;# create interpolation function
        range_acc_step_3 = range(1,stop=size_acc_step_3,length=200);
        interpolated_acc_step_3 = itp(range_acc_step_3);
        MAT_acc_step_3 = hcat(MAT_acc_step_3,interpolated_acc_step_3);

end

MAT_acc_step_1 = MAT_acc_step_1[:,2:end]
MAT_acc_step_2 = MAT_acc_step_2[:,2:end]
MAT_acc_step_3 = MAT_acc_step_3[:,2:end]

MAT_acc_step_1_right = MAT_acc_step_1[:,1:2:end]
MAT_acc_step_2_right = MAT_acc_step_2[:,1:2:end]
MAT_acc_step_3_right = MAT_acc_step_3[:,1:2:end]

MAT_acc_step_1_left = MAT_acc_step_1[:,2:2:end]
MAT_acc_step_2_left = MAT_acc_step_2[:,2:2:end]
MAT_acc_step_3_left = MAT_acc_step_3[:,2:2:end]


mean_of_size_acc_step_1_right = mean(MAT_acc_step_1_right,dims=2)
mean_of_size_acc_step_2_right = mean(MAT_acc_step_2_right,dims=2)
mean_of_size_acc_step_3_right = mean(MAT_acc_step_3_right,dims=2)

mean_of_size_acc_step_1_left = mean(MAT_acc_step_1_left,dims=2)
mean_of_size_acc_step_2_left = mean(MAT_acc_step_2_left,dims=2)
mean_of_size_acc_step_3_left = mean(MAT_acc_step_3_left,dims=2)



begin
    plt = plot(mean_of_size_acc_step_1_right)
	#plt = plot(acc)
	savefig("oussama_test_mean1_right.png")
end

begin
    plt = plot(mean_of_size_acc_step_2_right)
	#plt = plot(acc)
	savefig("oussama_test_mean2_right.png")
end

begin
    plt = plot(mean_of_size_acc_step_3_right)
	#plt = plot(acc)
	savefig("oussama_test_mean3_right.png")
end

begin
    plt = plot(mean_of_size_acc_step_1_left)
	#plt = plot(acc)
	savefig("oussama_test_mean1_left.png")
end

begin
    plt = plot(mean_of_size_acc_step_2_left)
	#plt = plot(acc)
	savefig("oussama_test_mean2_left.png")
end

begin
    plt = plot(mean_of_size_acc_step_3_left)
	#plt = plot(acc)
	savefig("oussama_test_mean3_left.png")
end

begin
    plt = plot(acc)
	#plt = plot(acc)
	savefig("oussama.png")
end