using Plots
using ActIO.File
r = recording("6025376C.hdr")
acc = r[accel,54500:55500]
acc_1 = acc[:,1]                                                                                                            
acc_2 = acc[:,2]
acc_3 = acc[:,3]

begin 
    plot(acc);
    savefig("plot_video.png")
end

using ActAnalysis.StepSlc
using DataFrames
using CSV
using ActAnalysis.Activity
using ActAnalysis.SpeedSvr
using ActAnalysis.Stepwave
using ActAnalysis.StepSlc
using Interpolations
using Statistics

stepslcdetector = StepSlcDetector()
slcsteps = stepslcdetector(acc)

DataFrame(slcsteps)

stepanalyser = StepAnalyser()
steps = stepanalyser(acc)
DataFrame(steps)

group_steps = steps[1:4,:]
size_group = 1:3;
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

#counter = counter+1;
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
sd_of_size_acc_step_1_right = std(MAT_acc_step_1_right,dims=2)
mean_of_size_acc_step_2_right = mean(MAT_acc_step_2_right,dims=2)
sd_of_size_acc_step_2_right = std(MAT_acc_step_2_right,dims=2)
mean_of_size_acc_step_3_right = mean(MAT_acc_step_3_right,dims=2)
sd_of_size_acc_step_3_right = std(MAT_acc_step_3_right,dims=2)


mean_of_size_acc_step_1_left = mean(MAT_acc_step_1_left,dims=2)
sd_of_size_acc_step_1_left = std(MAT_acc_step_1_left,dims=2)
mean_of_size_acc_step_2_left = mean(MAT_acc_step_2_left,dims=2)
sd_of_size_acc_step_2_left = std(MAT_acc_step_2_left,dims=2)
mean_of_size_acc_step_3_left = mean(MAT_acc_step_3_left,dims=2)
sd_of_size_acc_step_3_left = std(MAT_acc_step_3_left,dims=2)


normstep_variability_right=sqrt(sum( sd_of_size_acc_step_1_right.^ 2))/200 + sqrt(sum( sd_of_size_acc_step_2_right.^ 2))/200 + sqrt(sum( sd_of_size_acc_step_3_right.^ 2))/200 
normstep_variability_left=sqrt(sum( sd_of_size_acc_step_1_left.^ 2))/200 + sqrt(sum( sd_of_size_acc_step_2_left.^ 2))/200 + sqrt(sum( sd_of_size_acc_step_3_left.^ 2))/200 


diff1 = mean_of_size_acc_step_1_right - mean_of_size_acc_step_1_left;
diff2 = mean_of_size_acc_step_2_right - mean_of_size_acc_step_3_left;
diff3 = mean_of_size_acc_step_3_right - mean_of_size_acc_step_2_left;


for i = 1:200
    if diff1[i]<0
diff1[i]=-diff1[i]
end
end

for i = 1:200
    if diff2[i]<0
diff2[i]=-diff2[i]
end
end


for i = 1:200
    if diff3[i]<0
diff3[i]=-diff3[i]
end
end

normtep_asym = sum( diff1)/200 +sum( diff2)/200 + sum( diff3)/200 

acc_step_1_right_lower = mean_of_size_acc_step_1_right - sd_of_size_acc_step_1_right;
acc_step_1_right_upper = mean_of_size_acc_step_1_right + sd_of_size_acc_step_1_right;
acc_step_2_right_lower = mean_of_size_acc_step_2_right - sd_of_size_acc_step_2_right;
acc_step_2_right_upper = mean_of_size_acc_step_2_right + sd_of_size_acc_step_2_right;
acc_step_3_right_lower = mean_of_size_acc_step_3_right - sd_of_size_acc_step_3_right;
acc_step_3_right_upper = mean_of_size_acc_step_3_right + sd_of_size_acc_step_3_right;

acc_step_1_left_lower = mean_of_size_acc_step_1_left - sd_of_size_acc_step_1_left;
acc_step_1_left_upper = mean_of_size_acc_step_1_left + sd_of_size_acc_step_1_left;
acc_step_2_left_lower = mean_of_size_acc_step_2_left - sd_of_size_acc_step_2_left;
acc_step_2_left_upper = mean_of_size_acc_step_2_left + sd_of_size_acc_step_2_left;
acc_step_3_left_lower = mean_of_size_acc_step_3_left - sd_of_size_acc_step_3_left;
acc_step_3_left_upper = mean_of_size_acc_step_3_left + sd_of_size_acc_step_3_left;


title = plot(title = "Normstep", grid = false,showaxis=false,xaxis=nothing,yaxis=nothing, bottom_margin = -10Plots.px)
p1_L = plot(mean_of_size_acc_step_1_left,title = "Left foot",xlabel = "One Step samples", ylabel = "Acceleration",label = "x",lw = 1.5,color = :red)
plot!(acc_step_1_left_lower,title = "Left foot",xlabel = "One Step samples", ylabel = "Acceleration",label = "x-x.sd",lw = 0.1,color = :red)
plot!(acc_step_1_left_upper,title = "Left foot",xlabel = "One Step samples", ylabel = "Acceleration",label = "x+x.sd",lw = 0.1,color = :red)
plot!(mean_of_size_acc_step_2_left,title = "Left foot",xlabel = "One Step samples",label = "y",lw = 1.5,color = :blue)
plot!(acc_step_2_left_lower,title = "Left foot",xlabel = "One Step samples",label = "y-y.sd",lw = 0.1, color = :blue)
plot!(acc_step_2_left_upper,title = "Left foot",xlabel = "One Step samples",label = "y+y.sd",lw = 0.1, color = :blue)
plot!(mean_of_size_acc_step_3_left,title = "Left foot",xlabel = "One Step samples",label = "z",lw = 1.5,color = :green)
plot!(acc_step_3_left_lower,title = "Left foot",xlabel = "One Step samples",label = "z-z.sd",lw = 0.1,color = :green)
plot!(acc_step_3_left_upper,title = "Left foot",xlabel = "One Step samples",label = "z+z.sd",lw = 0.1,color = :green)
p1_R = plot(mean_of_size_acc_step_1_right,title = "Right foot",xlabel = "One Step samples", ylabel = "Acceleration",label = "x",lw = 1.5,color = :red)
plot!(acc_step_1_right_lower,title = "Right foot",xlabel = "One Step samples",label = "x-x.sd",lw = 0.1,color = :red)
plot!(acc_step_1_right_upper,title = "Right foot",xlabel = "One Step samples",label = "x+x.sd",lw = 0.1,color = :red)
plot!(mean_of_size_acc_step_2_right,title = "Right foot",xlabel = "One Step samples",label = "y",lw = 1.5,color = :blue)
plot!(acc_step_2_right_lower,title = "Right foot",xlabel = "One Step samples",label = "y-y.sd",lw = 0.1,color = :blue)
plot!(acc_step_2_right_upper,title = "Right foot",xlabel = "One Step samples",label = "y+y.sd",lw = 0.1,color = :blue)
plot!(mean_of_size_acc_step_3_right,title = "Right foot",xlabel = "One Step samples",label = "z",lw = 1.5,color = :green)
plot!(acc_step_3_right_lower,title = "Right foot",xlabel = "One Step samples",label = "z-z.sd",lw = 0.1,color = :green)
plot!(acc_step_3_right_upper,title = "Right foot",xlabel = "One Step samples",label = "z+z.sd",lw = 0.1,color = :green)
plot(title,p1_L, p1_R,layout = @layout([A{0.01h}; [B C]]))
savefig("aydin_left_right_acc_video.png")

