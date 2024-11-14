clc; clear all; close all;
%% Generating the Problem Table
filename = 'CDU.xlsx';%'CDU.xlsx';%'StreamData.xlsx';
streamData = readmatrix(filename);
Tmin = 5.8;
streamData(:,4) = streamData(:,1).*(streamData(:,2) - streamData(:,3));
streamData = streamData(:,1:4);

is_hot_stream = streamData(:,4) > 0;
hot_Tin_vec = streamData(is_hot_stream,2);
hot_Tout_vec = streamData(is_hot_stream,3);
cold_Tin_vec = streamData(~is_hot_stream,2);
cold_Tout_vec = streamData(~is_hot_stream,3);

Tin_vec = streamData(:,2);
Tout_vec = streamData(:,3);

hot_stream_mat = streamData(is_hot_stream,:);
cold_stream_mat = streamData(~is_hot_stream,:);

n_h = sum(is_hot_stream);
n_c = sum(~is_hot_stream);
n_tot = n_c + n_h;

prob_table = [[cold_Tin_vec+Tmin, cold_Tin_vec];[cold_Tout_vec+Tmin, cold_Tout_vec];...
              [hot_Tin_vec, hot_Tin_vec-Tmin];[hot_Tout_vec, hot_Tout_vec-Tmin]];
prob_table = unique(prob_table, 'rows');
prob_table = sortrows(prob_table, -1);

n_prob = length(prob_table);

% 1    2        3           4       5           6    7      8        9     10
% T_h | T_c | Avail heat | Casc | Reqd Heat | Casc | Net | cum Net | GCC | Adj Hot 
% Ignore last row --> For first interval check the first row and so on
prob_table(:,3:9) = 0;

for i = 1:n_h
    FC_i = hot_stream_mat(i,1);
    Tin_i = hot_stream_mat(i,2);
    Tout_i = hot_stream_mat(i,3);
    
    T_range_i = (prob_table(:,1) > Tout_i )&(prob_table(:,1) <= Tin_i );
    prob_table(1:end-1,3) = prob_table(1:end-1,3) + FC_i*(prob_table(1:end-1,1)-prob_table(2:end,1)).*T_range_i(1:end-1);
end
prob_table(2:end,4) = cumsum(prob_table(1:end-1,3));

for i = 1:n_c
    FC_i = cold_stream_mat(i,1);
    Tin_i = cold_stream_mat(i,3);
    Tout_i = cold_stream_mat(i,2);
    
    T_range_i = (prob_table(:,2) > Tout_i )&(prob_table(:,2) <= Tin_i );
    prob_table(1:end-1,5) = prob_table(1:end-1,5) + FC_i*(prob_table(1:end-1,2)-prob_table(2:end,2)).*T_range_i(1:end-1);
end
prob_table(2:end,6) = cumsum(prob_table(1:end-1,5));

prob_table(:,7) = prob_table(:,3) - prob_table(:,5);
prob_table(2:end,8) = cumsum(prob_table(1:end-1,7));

if min(prob_table(:,8)) < 0
    h_utility = min(prob_table(:,8));
    prob_table(:,9) = prob_table(:,8) - h_utility;
    prob_table(:,10) = prob_table(:,4) - h_utility;
else
    h_utility = 0;
end

[c_h,ia_h,~] = unique(prob_table(:,10), 'rows');
[c_c,ia_c,~] = unique(prob_table(:,6), 'rows');

f1 = figure;
plot(c_c,prob_table(ia_c,2),'-*','linewidth',1.5)
hold on
plot(c_h,prob_table(ia_h,1),'-*','linewidth',1.5)
xlabel('Cascaded Heat (kW)')
ylabel('Temperature (K)')
title('TQ Curve')

f2 = figure;
plot(prob_table(:,9),prob_table(:,1),'-*','linewidth',1.5)
xlabel('Cascaded Heat (kW)')
ylabel('Temperature (K)')
title('Grand Composite Curve')

%% Miniumm HEx
% Obtaining initial solution
pinch_temp_h = prob_table(prob_table(:,9)==0,1);

is_above_pinch = (streamData(:,2) > pinch_temp_h) | (streamData(:,3) > pinch_temp_h-Tmin) ;
is_below_pinch = (streamData(:,2) < pinch_temp_h-Tmin) | (streamData(:,3) < pinch_temp_h);

% Above pinch streams
% Heat in each hot stream and hot utility at the end

%% Above Pinch
heat_available_above = [streamData(is_above_pinch & is_hot_stream,1).*(streamData(is_above_pinch & is_hot_stream,2)-max(streamData(is_above_pinch & is_hot_stream,3),pinch_temp_h)); prob_table(1,9)];  
heat_available_above_info = streamData(is_above_pinch & is_hot_stream,:);
heat_available_above_info(end+1,4) = prob_table(1,9);
heat_required_above = [(-1).*streamData(is_above_pinch & ~is_hot_stream,1).*(max(streamData(is_above_pinch & ~is_hot_stream,2),pinch_temp_h-Tmin)-streamData(is_above_pinch & ~is_hot_stream,3))];  
heat_required_above_info = streamData(is_above_pinch & ~is_hot_stream,:);

num_hot_above = length(heat_available_above);
num_cold_above = length(heat_required_above);

% Initialize adjacency matrix for above pinch
adj_matrix_above = zeros(num_hot_above, num_cold_above);

% Populate adjacency matrix with heat transfer estimates for above pinch streams
for i = 1:num_hot_above
    for j = 1:num_cold_above
        adj_matrix_above(i, j) = min(heat_available_above(i), heat_required_above(j));
        heat_available_above(i) = heat_available_above(i) - adj_matrix_above(i, j);
        heat_required_above(j) = heat_required_above(j) - adj_matrix_above(i, j);
    end
end

% Display initial adjacency matrix for above pinch
disp('Initial Adjacency Matrix (Above Pinch):');
disp(adj_matrix_above);

% Below Pinch
% Define heat availability for hot streams below the pinch (hot utility) and heat requirements for cold streams below pinch
heat_available_below = [streamData(is_below_pinch & is_hot_stream,1).*(min(streamData(is_below_pinch & is_hot_stream,2),pinch_temp_h)-streamData(is_below_pinch & is_hot_stream,3))];  
heat_available_below_info = streamData(is_below_pinch & is_hot_stream,:);
heat_required_below = [(-1).*streamData(is_below_pinch & ~is_hot_stream,1).*(streamData(is_below_pinch & ~is_hot_stream,2)-min(streamData(~is_hot_stream & is_below_pinch,3),pinch_temp_h-Tmin)); prob_table(end,9)]; 
heat_required_below_info = streamData(is_below_pinch & ~is_hot_stream,:);
heat_required_below_info(end+1,4) = prob_table(end,9);

num_hot_below = length(heat_available_below);
num_cold_below = length(heat_required_below);

% Initialize adjacency matrix for below pinch
adj_matrix_below = zeros(num_hot_below, num_cold_below);

% Populate adjacency matrix with heat transfer estimates for below pinch streams
for i = 1:num_hot_below
    for j = 1:num_cold_below
        adj_matrix_below(i, j) = min(heat_available_below(i), heat_required_below(j));
        heat_available_below(i) = heat_available_below(i) - adj_matrix_below(i, j);
        heat_required_below(j) = heat_required_below(j) - adj_matrix_below(i, j);
    end
end

% Display initial adjacency matrix for below pinch
disp('Initial Adjacency Matrix (Below Pinch):');
disp(adj_matrix_below);



