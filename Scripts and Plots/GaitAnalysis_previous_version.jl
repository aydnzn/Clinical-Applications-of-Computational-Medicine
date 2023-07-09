## Aydin Uzun & Oussama Skhiri 
# Using these packages
using ActIO.File
using ActAnalysis.Stepwave
using DataFrames
using Dates
using Plots
using Statistics
using Interpolations
using LinearAlgebra
using JLD2

# individual number
idx = 1
# normstep initialization
Normstep_acc_x_L_prev = zeros(1, 100);
Normstep_acc_y_L_prev = zeros(1, 100);
Normstep_acc_z_L_prev = zeros(1, 100);
Normstep_acc_x_R_prev = zeros(1, 100);
Normstep_acc_y_R_prev = zeros(1, 100);
Normstep_acc_z_R_prev = zeros(1, 100);
# step analyser function from actibelt package
sa = StepAnalyser(300, 0.8, (6, 6, 6, 21), (0.01, 0.8), (0.02, 0.65), (0.15, 0.7), (0.75, 0.95), (100, 0.2), 0.1, 290);
# get our own data
r_aydin_3 = recording("60251AD6.HDR");
r_aydin_2 = recording("601B568A.HDR");
r_oussama = recording("oussama.HDR");
# get the data from actibelt database
@load "walking_data_extracts.jld2" acceleration_data_DW

# concatanate them
array_acc = vcat([r_aydin_2[accel],r_aydin_3[accel],r_oussama[accel]], acceleration_data_DW)
L2_Norm_left = zeros(3, 26)
L2_Norm_right = zeros(3, 26)

# for plot enhancement
gr(size=(1200, 1000), xtickfontsize=13, ytickfontsize=13, xguidefontsize=16, yguidefontsize=16, legendfontsize=10, dpi=100, grid=(:y, :gray, :solid, 1, 0.4));

for acc in array_acc
    # get x,y,z accelerations
    acc_z = acc[:,1];
    acc_y = acc[:,2];
    acc_x = acc[:,3];
    # get steps
    steps = (sa)(acc, frequency_extractor=FrequencyExtractor(), level=5);
    # get index of left and right foot in the array
    idx_right_foot = findall(x -> x == "right", steps.side);
    idx_left_foot = findall(x -> x == "left", steps.side);
    # create coherent step groups ...
    # eliminate steps whose duration not in [0.3 1]
    steps_left = steps[idx_left_foot,:];
    steps_right = steps[idx_right_foot,:];
    idx_left_foot = findall(x -> (x < 1 && x > 0.3), steps_left.duration);
    idx_right_foot = findall(x -> (x < 1 && x > 0.3), steps_right.duration);
    steps_left_filtered = steps_left[idx_left_foot,:];
    steps_right_filtered = steps_right[idx_right_foot,:];
    # get the interpolation size which corresponds to the largest duration converted in samples
    # interp = maximum([steps_left_filtered.duration;steps_right_filtered.duration]);
    interp = 1;
    interp_samples = Int64(round(interp * 100));    
    # Normstep for the left foot.
    MAT_acc_z_step_L = zeros(size(steps_left_filtered, 1), interp_samples);
    MAT_acc_y_step_L = zeros(size(steps_left_filtered, 1), interp_samples);
    MAT_acc_x_step_L = zeros(size(steps_left_filtered, 1), interp_samples);
    for i in 1:size(steps_left_filtered, 1)
        # get time of begining of the step + its duration
        time = steps_left_filtered[i,1];
        dur = steps_left_filtered[i,4];
        # get index of start and end timestep in the acc data
        start_step_idx = datetime_to_index(time, start(acc), 10) - 1;
        end_step_idx = datetime_to_index(time + Dates.Millisecond(round(dur * 1000)), start(acc), 10) - 1;
        acc_z_in_step = acc_z[start_step_idx:end_step_idx];
        acc_y_in_step = acc_y[start_step_idx:end_step_idx];
        acc_x_in_step = acc_x[start_step_idx:end_step_idx];
        # interpolate the array to get array with length equals to inter_samples
        itp = LinearInterpolation(1:size(acc_z_in_step, 1), acc_z_in_step) ;
        range_acc_step = range(1, stop=size(acc_z_in_step, 1), length=interp_samples);
        interpolated_acc_z_step = itp(range_acc_step);    
        # interpolate the array to get array with length equals to inter_samples
        itp = LinearInterpolation(1:size(acc_y_in_step, 1), acc_y_in_step) ;
        range_acc_step = range(1, stop=size(acc_y_in_step, 1), length=interp_samples);
        interpolated_acc_y_step = itp(range_acc_step);    
        # interpolate the array to get array with length equals to inter_samples
        itp = LinearInterpolation(1:size(acc_x_in_step, 1), acc_x_in_step) ;
        range_acc_step = range(1, stop=size(acc_x_in_step, 1), length=interp_samples);
        interpolated_acc_x_step = itp(range_acc_step);    
        # save in matrix
        MAT_acc_z_step_L[i,:] = interpolated_acc_z_step;
        MAT_acc_y_step_L[i,:] = interpolated_acc_y_step;
        MAT_acc_x_step_L[i,:] = interpolated_acc_x_step;
    end
    Normstep_acc_z_L = mean(MAT_acc_z_step_L, dims=1);
    Normstep_acc_y_L = mean(MAT_acc_y_step_L, dims=1);
    Normstep_acc_x_L = mean(MAT_acc_x_step_L, dims=1);
    
    # Normstep for the right foot
    MAT_acc_z_step_R = zeros(size(steps_right_filtered, 1), interp_samples);
    MAT_acc_y_step_R = zeros(size(steps_right_filtered, 1), interp_samples);
    MAT_acc_x_step_R = zeros(size(steps_right_filtered, 1), interp_samples);    
    for i in 1:size(steps_right_filtered, 1)
        # get time of begining of the step + its duration
        time_R = steps_right_filtered[i,1];
        dur_R = steps_right_filtered[i,4];
        # get index of start and end timestep in the acc data
        start_step_idx_R = datetime_to_index(time_R, start(acc), 10) - 1;
        end_step_idx_R = datetime_to_index(time_R + Dates.Millisecond(round(dur_R * 1000)), start(acc), 10) - 1;
        acc_z_in_step_R = acc_z[start_step_idx_R:end_step_idx_R];
        acc_y_in_step_R = acc_y[start_step_idx_R:end_step_idx_R];
        acc_x_in_step_R = acc_x[start_step_idx_R:end_step_idx_R];
        # interpolate the array to get array with length equals to inter_samples
        itp_R = LinearInterpolation(1:size(acc_z_in_step_R, 1), acc_z_in_step_R) ;
        range_acc_step_R = range(1, stop=size(acc_z_in_step_R, 1), length=interp_samples);
        interpolated_acc_z_step_R = itp_R(range_acc_step_R);    
        # interpolate the array to get array with length equals to inter_samples
        itp_R = LinearInterpolation(1:size(acc_y_in_step_R, 1), acc_y_in_step_R) ;
        range_acc_step_R = range(1, stop=size(acc_y_in_step_R, 1), length=interp_samples);
        interpolated_acc_y_step_R = itp_R(range_acc_step_R);    
        # interpolate the array to get array with length equals to inter_samples
        itp_R = LinearInterpolation(1:size(acc_x_in_step_R, 1), acc_x_in_step_R) ;
        range_acc_step_R = range(1, stop=size(acc_x_in_step_R, 1), length=interp_samples);
        interpolated_acc_x_step_R = itp_R(range_acc_step_R);    
        # save in matrix
        MAT_acc_z_step_R[i,:] = interpolated_acc_z_step_R;
        MAT_acc_y_step_R[i,:] = interpolated_acc_y_step_R;
        MAT_acc_x_step_R[i,:] = interpolated_acc_x_step_R;
    end
    Normstep_acc_z_R = mean(MAT_acc_z_step_R, dims=1);
    Normstep_acc_y_R = mean(MAT_acc_y_step_R, dims=1);
    Normstep_acc_x_R = mean(MAT_acc_x_step_R, dims=1);

    # calculate the standard deviations
    sd_Normstep_acc_z_R = std(MAT_acc_z_step_R, dims=1);
    sd_Normstep_acc_y_R = std(MAT_acc_y_step_R, dims=1);
    sd_Normstep_acc_x_R = std(MAT_acc_x_step_R, dims=1);

    sd_Normstep_acc_z_L = std(MAT_acc_z_step_L, dims=1);
    sd_Normstep_acc_y_L = std(MAT_acc_y_step_L, dims=1);
    sd_Normstep_acc_x_L = std(MAT_acc_x_step_L, dims=1);

    acc_z_R_lower = Normstep_acc_z_R - sd_Normstep_acc_z_R;
    acc_z_R_upper = Normstep_acc_z_R + sd_Normstep_acc_z_R;
    acc_y_R_lower = Normstep_acc_y_R - sd_Normstep_acc_y_R;
    acc_y_R_upper = Normstep_acc_y_R + sd_Normstep_acc_y_R;
    acc_x_R_lower = Normstep_acc_x_R - sd_Normstep_acc_x_R;
    acc_x_R_upper = Normstep_acc_x_R + sd_Normstep_acc_x_R;

    acc_z_L_lower = Normstep_acc_z_L - sd_Normstep_acc_z_L;
    acc_z_L_upper = Normstep_acc_z_L + sd_Normstep_acc_z_L;
    acc_y_L_lower = Normstep_acc_y_L - sd_Normstep_acc_y_L;
    acc_y_L_upper = Normstep_acc_y_L + sd_Normstep_acc_y_L;
    acc_x_L_lower = Normstep_acc_x_L - sd_Normstep_acc_x_L;
    acc_x_L_upper = Normstep_acc_x_L + sd_Normstep_acc_x_L;

    # normstep plot
    title = plot(title="Normstep", grid=false, showaxis=false, xaxis=nothing, yaxis=nothing, bottom_margin=-10Plots.px);
    p1_L = plot(Normstep_acc_x_L', title="Left foot", xlabel="One Step samples", ylabel="Acceleration", label="x", lw=1.5, color=:red, ylims=(-2, 1.5));
    plot!(acc_x_L_lower', fillrange=acc_x_L_upper',fillalpha=0.2,c=:red,title="Left foot",xlabel="One Step samples", ylabel="Acceleration",label="confidence band",lw=0.5,color=:red );
    plot!(Normstep_acc_y_L',title="Left foot",xlabel="One Step samples", ylabel="Acceleration",label="y",lw=1.5,color=:blue,  ylims=(-2, 1.5) );
    plot!(acc_y_L_lower', fillrange=acc_y_L_upper',fillalpha=0.2,c=:blue,title="Left foot",xlabel="One Step samples", ylabel="Acceleration",label="confidence band",lw=0.5,color=:blue );
    plot!(Normstep_acc_z_L',title="Left foot",xlabel="One Step samples", ylabel="Acceleration",label="z",lw=1.5,color=:green , ylims=(-2, 1.5));
    plot!(acc_z_L_lower', fillrange=acc_z_L_upper',fillalpha=0.2,c=:green,title="Left foot",xlabel="One Step samples", ylabel="Acceleration",label="confidence band",lw=0.5,color=:green );
    p1_R = plot(Normstep_acc_x_R', title="Right foot", xlabel="One Step samples", ylabel="Acceleration", label="x", lw=1.5, color=:red, ylims=(-2, 1.5));
    plot!(acc_x_R_lower', fillrange=acc_x_R_upper',fillalpha=0.2,c=:red,title="Right foot",xlabel="One Step samples", ylabel="Acceleration",label="confidence band",lw=0.5,color=:red );
    plot!(Normstep_acc_y_R',title="Right foot",xlabel="One Step samples", ylabel="Acceleration",label="y",lw=1.5,color=:blue,ylims=(-2, 1.5));
    plot!(acc_y_R_lower', fillrange=acc_y_R_upper',fillalpha=0.2,c=:blue,title="Right foot",xlabel="One Step samples", ylabel="Acceleration",label="confidence band",lw=0.5,color=:blue );
    plot!(Normstep_acc_z_R',title="Right foot",xlabel="One Step samples", ylabel="Acceleration",label="z",lw=1.5,color=:green,ylims=(-2, 1.5));
    plot!(acc_z_R_lower', fillrange=acc_z_R_upper',fillalpha=0.2,c=:green,title="Right foot",xlabel="One Step samples", ylabel="Acceleration",label="confidence band",lw=0.5,color=:green );
    plot(title,p1_L, p1_R,layout=@layout([A{0.01h}; [B C]]));
    savefig("NormStep" * string(idx) * ".png");
    # L2 Norm of the difference        
    diff_x_L = Normstep_acc_x_L - Normstep_acc_x_L_prev;
    diff_y_L = Normstep_acc_y_L - Normstep_acc_y_L_prev;
    diff_z_L = Normstep_acc_z_L - Normstep_acc_z_L_prev;
    diff_x_R = Normstep_acc_x_R - Normstep_acc_x_R_prev;
    diff_y_R = Normstep_acc_y_R - Normstep_acc_y_R_prev;
    diff_z_R = Normstep_acc_z_R - Normstep_acc_z_R_prev;
    # save the first normstep for oussama to compare it with another excerpt of oussama and also other Persons
    if (idx == 1)
        Normstep_acc_x_L_prev = Normstep_acc_x_L;
        Normstep_acc_y_L_prev = Normstep_acc_y_L;
        Normstep_acc_z_L_prev = Normstep_acc_z_L;
        Normstep_acc_x_R_prev = Normstep_acc_x_R;
        Normstep_acc_y_R_prev = Normstep_acc_y_R;
        Normstep_acc_z_R_prev = Normstep_acc_z_R;
    end
    L2_Norm_left[:,idx] = [norm(diff_x_L, 2);norm(diff_y_L, 2);norm(diff_z_L, 2)];
    L2_Norm_right[:,idx] = [norm(diff_x_R, 2);norm(diff_y_R, 2);norm(diff_z_R, 2)];
    idx += 1;
end 
# Plot L2 Norm for left and right steps for the three acceleration parts
title = plot(title="Comparison of Normstep excerpts between different individuals", grid=false, showaxis=false, xaxis=nothing, yaxis=nothing, bottom_margin=-10Plots.px)
p1_L = plot(L2_Norm_left[1,2:end], title="Left foot", seriestype=:scatter, xaxis=("Comparisons", (0, 26), 1:25), xrotation=70, ylabel="L2 Norm", label="x", legend=:topright, ylims=(0, 10) )
plot!(L2_Norm_left[2,2:end],title="Left foot",xaxis=("Comparisons", (0, 26), 1:25),seriestype=:scatter,xrotation=70,label="y",legend=:topright)# ,size=(850,800)
plot!(L2_Norm_left[3,2:end],title="Left foot",xaxis=("Comparisons", (0, 26), 1:25),xrotation=70,seriestype=:scatter,label="z",legend=:topright)
p1_R = plot(L2_Norm_right[1,2:end], title="Right foot", xaxis=("Comparisons", (0, 26), 1:25), xrotation=70, seriestype=:scatter, ylabel="L2 Norm", label="x", legend=:topright, ylims=(0,10) )
plot!(L2_Norm_right[2,2:end],title="Right foot",xaxis=("Comparisons", (0, 26), 1:25),xrotation=70,seriestype=:scatter,label="y",legend=:topright)
plot!(L2_Norm_right[3,2:end],title="Right foot",xaxis=("Comparisons", (0, 26), 1:25),xrotation=70,seriestype=:scatter,label="z",legend=:topright)
plot(title,p1_L, p1_R,layout=@layout([A{0.01h}; [B C]]))
savefig("L2Norm_Diff.png")




