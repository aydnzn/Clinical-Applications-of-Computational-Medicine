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
using StatsBase
# step analyser function from actibelt package
sa = StepAnalyser();

# get our own data
r_aydin_3 = recording("60251AD6.HDR");
r_aydin_2 = recording("601B568A.HDR");
r_oussama = recording("oussama.HDR");
# get the data from actibelt database
@load "walking_data_extracts.jld2" acceleration_data_DW
@load "walking_data_extracts_part2.jld2" acceleration_data_DW_part2

# concatanate them
array_acc = vcat([r_aydin_2[accel],r_aydin_3[accel],r_oussama[accel]], acceleration_data_DW, acceleration_data_DW_part2)

# initialize normstep matrices
MAT_norsmtep_z_R = zeros(100, 1);
MAT_norsmtep_y_R = zeros(100, 1);
MAT_norsmtep_x_R = zeros(100, 1);
MAT_norsmtep_z_L = zeros(100, 1);
MAT_norsmtep_y_L = zeros(100, 1);
MAT_norsmtep_x_L = zeros(100, 1);

# for plot enhancement
gr(size=(1200, 1000), xtickfontsize=13, ytickfontsize=13, xguidefontsize=16, yguidefontsize=16, legendfontsize=10, dpi=100, grid=(:y, :gray, :solid, 1, 0.4));

idx =1; # for plotting
for acc in array_acc
    # get x,y,z accelerations
    acc_z = acc[:,1];
    acc_y = acc[:,2];
    acc_x = acc[:,3];
    # get steps
    steps = sa(acc);
    # get index of left and right foot in the array
    idx_right_foot = findall(x -> x == "right", steps.side);
    idx_left_foot = findall(x -> x == "left", steps.side);
    # create coherent step groups ...
    # eliminate steps whose duration not in [0.3 1]
    steps_left = steps[idx_left_foot,:];
    steps_right = steps[idx_right_foot,:];
    idx_left_foot = findall(x -> (x < 1 && x > 0.3), steps_left.duration);
    idx_right_foot = findall(x -> (x < 1 && x > 0.3), steps_right.duration);
    steps_left_filtered_before_bout = steps_left[idx_left_foot,:];
    steps_right_filtered_before_bout = steps_right[idx_right_foot,:];
    # eliminate steps with speeds not in the range [0.55, 1.5]
    idx_left_foot_speed = findall(x -> (x > 0.55 && x < 1.5), steps_left_filtered_before_bout.speed);
    idx_right_foot_speed = findall(x -> (x > 0.55 && x < 1.5), steps_right_filtered_before_bout.speed);
    steps_left_filtered_before_bout_speed = steps_left_filtered_before_bout[idx_left_foot_speed,:];
    steps_right_filtered_before_bout_speed = steps_right_filtered_before_bout[idx_right_foot_speed,:];
    # count the uninterrupted steps in the filtered step data 
    countmap_bout = countmap(steps_left_filtered_before_bout_speed.bout);
    # find maximum of uninterrupted steps
    bout_and_occurrence = findmax(countmap_bout);
    max_occurring_bout = bout_and_occurrence[2];
    occurrence = bout_and_occurrence[1];

    # if filtered uninterrupted step group larger than 40...
    if occurrence > 40 
        idx_left_foot_group = findall(x -> x == max_occurring_bout, steps_left_filtered_before_bout_speed.bout);
        idx_right_foot_group = findall(x -> x == max_occurring_bout, steps_right_filtered_before_bout_speed.bout);
        steps_left_filtered = steps_left_filtered_before_bout_speed[idx_left_foot_group,:];
        steps_right_filtered = steps_right_filtered_before_bout_speed[idx_right_foot_group,:];
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
        # calculate the normstep
        Normstep_acc_z_L = mean(MAT_acc_z_step_L, dims=1);
        Normstep_acc_y_L = mean(MAT_acc_y_step_L, dims=1);
        Normstep_acc_x_L = mean(MAT_acc_x_step_L, dims=1);

        # concatanate the normsteps in a matrix
        MAT_norsmtep_z_L = hcat(MAT_norsmtep_z_L, Normstep_acc_z_L');
        MAT_norsmtep_y_L = hcat(MAT_norsmtep_y_L, Normstep_acc_y_L');
        MAT_norsmtep_x_L = hcat(MAT_norsmtep_x_L, Normstep_acc_x_L');

        ###### do the same for right foot
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
        # calculate the normstep
        Normstep_acc_z_R = mean(MAT_acc_z_step_R, dims=1);
        Normstep_acc_y_R = mean(MAT_acc_y_step_R, dims=1);
        Normstep_acc_x_R = mean(MAT_acc_x_step_R, dims=1);

        # concatanate the normsteps in a matrix
        MAT_norsmtep_z_R = hcat(MAT_norsmtep_z_R, Normstep_acc_z_R');
        MAT_norsmtep_y_R = hcat(MAT_norsmtep_y_R, Normstep_acc_y_R');
        MAT_norsmtep_x_R = hcat(MAT_norsmtep_x_R, Normstep_acc_x_R');

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

        # normstep plot for each individual
        title = plot(title="Normstep", grid=false, showaxis=false, xaxis=nothing, yaxis=nothing, bottom_margin=-10Plots.px);
        p1_L = plot(Normstep_acc_x_L', title="Left foot", xlabel="One Step samples", ylabel="Acceleration", label="x", lw=1.5, color=:red, ylims=(-2, 1.5));
        plot!(acc_x_L_lower', fillrange=acc_x_L_upper', fillalpha=0.2, c=:red, title="Left foot", xlabel="One Step samples", ylabel="Acceleration", label="confidence band", lw=0.5, color=:red);
        plot!(Normstep_acc_y_L', title="Left foot", xlabel="One Step samples", ylabel="Acceleration", label="y", lw=1.5, color=:blue,  ylims=(-2, 1.5));
        plot!(acc_y_L_lower', fillrange=acc_y_L_upper', fillalpha=0.2, c=:blue, title="Left foot", xlabel="One Step samples", ylabel="Acceleration", label="confidence band", lw=0.5, color=:blue);
        plot!(Normstep_acc_z_L', title="Left foot", xlabel="One Step samples", ylabel="Acceleration", label="z", lw=1.5, color=:green, ylims=(-2, 1.5));
        plot!(acc_z_L_lower', fillrange=acc_z_L_upper', fillalpha=0.2, c=:green, title="Left foot", xlabel="One Step samples", ylabel="Acceleration", label="confidence band", lw=0.5, color=:green);
        p1_R = plot(Normstep_acc_x_R', title="Right foot", xlabel="One Step samples", ylabel="Acceleration", label="x", lw=1.5, color=:red, ylims=(-2, 1.5));
        plot!(acc_x_R_lower', fillrange=acc_x_R_upper', fillalpha=0.2, c=:red, title="Right foot", xlabel="One Step samples", ylabel="Acceleration", label="confidence band", lw=0.5, color=:red);
        plot!(Normstep_acc_y_R', title="Right foot", xlabel="One Step samples", ylabel="Acceleration", label="y", lw=1.5, color=:blue, ylims=(-2, 1.5));
        plot!(acc_y_R_lower', fillrange=acc_y_R_upper', fillalpha=0.2, c=:blue, title="Right foot", xlabel="One Step samples", ylabel="Acceleration", label="confidence band", lw=0.5, color=:blue);
        plot!(Normstep_acc_z_R', title="Right foot", xlabel="One Step samples", ylabel="Acceleration", label="z", lw=1.5, color=:green, ylims=(-2, 1.5));
        plot!(acc_z_R_lower', fillrange=acc_z_R_upper', fillalpha=0.2, c=:green, title="Right foot", xlabel="One Step samples", ylabel="Acceleration", label="confidence band", lw=0.5, color=:green);
        plot(title, p1_L, p1_R, layout=@layout([A{0.01h}; [B C]]));
        savefig("NormStep" * string(idx) * ".png");
        idx += 1;
    end
end 

################ exclude first normstep because it is 0 vector by definition
Normstep_MAT_x_R  = MAT_norsmtep_x_R[:,2:end];
Normstep_MAT_y_R  = MAT_norsmtep_y_R[:,2:end];
Normstep_MAT_z_R  = MAT_norsmtep_z_R[:,2:end];
Normstep_MAT_x_L  = MAT_norsmtep_x_L[:,2:end]; 
Normstep_MAT_y_L  = MAT_norsmtep_y_L[:,2:end]; 
Normstep_MAT_z_L  = MAT_norsmtep_z_L[:,2:end]; 

#################### DONT RUN THIS PART - THIS PART IS MANUAL MODIFICATION TO EXCLUDE INDIVIDUAL NUMBER 11, 5 and 3 
# becasue they result in strong differences
Normstep_MAT_x_R = Normstep_MAT_x_R[:, 1:end .!= 11];
Normstep_MAT_x_R = Normstep_MAT_x_R[:, 1:end .!= 5];
Normstep_MAT_x_R = Normstep_MAT_x_R[:, 1:end .!= 3];

Normstep_MAT_y_R = Normstep_MAT_y_R[:, 1:end .!= 11];
Normstep_MAT_y_R = Normstep_MAT_y_R[:, 1:end .!= 5];
Normstep_MAT_y_R = Normstep_MAT_y_R[:, 1:end .!= 3];

Normstep_MAT_z_R = Normstep_MAT_z_R[:, 1:end .!= 11];
Normstep_MAT_z_R = Normstep_MAT_z_R[:, 1:end .!= 5];
Normstep_MAT_z_R = Normstep_MAT_z_R[:, 1:end .!= 3];

Normstep_MAT_x_L = Normstep_MAT_x_L[:, 1:end .!= 11];
Normstep_MAT_x_L = Normstep_MAT_x_L[:, 1:end .!= 5];
Normstep_MAT_x_L = Normstep_MAT_x_L[:, 1:end .!= 3];

Normstep_MAT_y_L = Normstep_MAT_y_L[:, 1:end .!= 11];
Normstep_MAT_y_L = Normstep_MAT_y_L[:, 1:end .!= 5];
Normstep_MAT_y_L = Normstep_MAT_y_L[:, 1:end .!= 3];

Normstep_MAT_z_L = Normstep_MAT_z_L[:, 1:end .!= 11];
Normstep_MAT_z_L = Normstep_MAT_z_L[:, 1:end .!= 5];
Normstep_MAT_z_L = Normstep_MAT_z_L[:, 1:end .!= 3];
####################################################


################## 36 is the individual number, if you have a different individual number in the end you should change it manually
# calculate L2 distance between each individual and store them in a matrice of size [individual number x individual number] 
compare_mat_x_R =  zeros(36, 36);
compare_mat_x_R_dummy = zeros(1, 36);
for i = 1:36
    ref_normstep = Normstep_MAT_x_R[:,i];
    to_be_subtracted = repeat(ref_normstep, 1, 36);
    distance_mat = Normstep_MAT_x_R - to_be_subtracted;

    for j = 1:36
        compare_mat_x_R_dummy[j] = norm(distance_mat[:,j]);
    end
    compare_mat_x_R[:,i] = compare_mat_x_R_dummy;
end

compare_mat_y_R =  zeros(36, 36);
compare_mat_y_R_dummy = zeros(1, 36);
for i = 1:36
    ref_normstep = Normstep_MAT_y_R[:,i];
    to_be_subtracted = repeat(ref_normstep, 1, 36);
    distance_mat = Normstep_MAT_y_R - to_be_subtracted;

    for j = 1:36
        compare_mat_y_R_dummy[j] = norm(distance_mat[:,j]);
    end
    compare_mat_y_R[:,i] = compare_mat_y_R_dummy;
end

compare_mat_z_R =  zeros(36, 36);
compare_mat_z_R_dummy = zeros(1, 36);
for i = 1:36
    ref_normstep = Normstep_MAT_z_R[:,i];
    to_be_subtracted = repeat(ref_normstep, 1, 36);
    distance_mat = Normstep_MAT_z_R - to_be_subtracted;

    for j = 1:36
        compare_mat_z_R_dummy[j] = norm(distance_mat[:,j]);
    end
    compare_mat_z_R[:,i] = compare_mat_z_R_dummy;
end

compare_mat_x_L =  zeros(36, 36);
compare_mat_x_L_dummy = zeros(1, 36);
for i = 1:36
    ref_normstep = Normstep_MAT_x_L[:,i];
    to_be_subtracted = repeat(ref_normstep, 1, 36);
    distance_mat = Normstep_MAT_x_L - to_be_subtracted;

    for j = 1:36
        compare_mat_x_L_dummy[j] = norm(distance_mat[:,j]);
    end
    compare_mat_x_L[:,i] = compare_mat_x_L_dummy;
end

compare_mat_y_L =  zeros(36, 36);
compare_mat_y_L_dummy = zeros(1, 36);
for i = 1:36
    ref_normstep = Normstep_MAT_y_L[:,i];
    to_be_subtracted = repeat(ref_normstep, 1, 36);
    distance_mat = Normstep_MAT_y_L - to_be_subtracted;

    for j = 1:36
        compare_mat_y_L_dummy[j] = norm(distance_mat[:,j]);
    end
    compare_mat_y_L[:,i] = compare_mat_y_L_dummy;
end

compare_mat_z_L =  zeros(36, 36);
compare_mat_z_L_dummy = zeros(1, 36);
for i = 1:36
    ref_normstep = Normstep_MAT_z_L[:,i];
    to_be_subtracted = repeat(ref_normstep, 1, 36);
    distance_mat = Normstep_MAT_z_L - to_be_subtracted;

    for j = 1:36
        compare_mat_z_L_dummy[j] = norm(distance_mat[:,j]);
    end
    compare_mat_z_L[:,i] = compare_mat_z_L_dummy;
end

# Plot the L2 distances as a heatmap
gr()
heatmap(1:36, 1:36, compare_mat_x_R, c=cgrad([:white,:black]), xlabel="Individuals", ylabel="Individuals", title="Normstep Right X Acceleration ");
savefig("normstep_R_x_acc_comparison.png")

gr()
heatmap(1:36, 1:36, compare_mat_y_R, c=cgrad([:white,:black]), xlabel="Individuals", ylabel="Individuals", title="Normstep Right Y Acceleration ");
savefig("normstep_R_y_acc_comparison.png")

gr()
heatmap(1:36, 1:36, compare_mat_z_R, c=cgrad([:white,:black]), xlabel="Individuals", ylabel="Individuals", title="Normstep Right Z Acceleration ");
savefig("normstep_R_z_acc_comparison.png")

gr()
heatmap(1:36, 1:36, compare_mat_x_L, c=cgrad([:white,:black]), xlabel="Individuals", ylabel="Individuals", title="Normstep Left X Acceleration ");
savefig("normstep_L_x_acc_comparison.png")

gr()
heatmap(1:36, 1:36, compare_mat_y_L, c=cgrad([:white,:black]), xlabel="Individuals", ylabel="Individuals", title="Normstep Left Y Acceleration ");
savefig("normstep_L_y_acc_comparison.png")

gr()
heatmap(1:36, 1:36, compare_mat_z_L, c=cgrad([:white,:black]), xlabel="Individuals", ylabel="Individuals", title="Normstep Left Z Acceleration ");
savefig("normstep_L_z_acc_comparison.png")
