load('/nas/longleaf/home/mrcole/ScratchingTheSurface/Code/exampleData/selected_subs_1000_totalcog.mat')

output_dir = "/nas/longleaf/home/mrcole/Encore/results/abcd_subject_lists/";

rng(544648654)

%% subjects for biological sex prediction
sex_0 = find(demographx(:,2) == 0);
sex_1 = find(demographx(:,2) == 1);

% shuffle the indices
sex_0 = sex_0(randperm(length(sex_0),length(sex_0)));
sex_1 = sex_1(randperm(length(sex_1),length(sex_1)));

% save the populations of interest
for N = [50,100,150,200,400]
    sub_idx = [sex_0(1:(N/2));sex_1(1:(N/2))];
    save(sprintf("%s/sex_N%0.3i.mat",output_dir,N), "sub_idx", "-v7.3");
end

%% subjects for cognitive analysis (both sexes)
% no cognitive scores for first participlant, so remove them
demographx(1,2) = 3;

for N = [50,100,150,200,400]
    sub_idx = [find(demographx(:,2) ~= 3, N/2, 'first')
               find(demographx(:,2) ~= 3, N/2, 'last')];

    save(sprintf("%s/cog_N%0.3i.mat",output_dir,N), "sub_idx", "-v7.3");
end

%% subjects for cognitive analysis (single sex)
for N = [50,100,150,200,400]
    sub_idx = [find(demographx(:,2) == 0, N/2, 'first')
               find(demographx(:,2) == 0, N/2, 'last')];

    save(sprintf("%s/cog_singlesex_N%0.3i.mat",output_dir,N), "sub_idx", "-v7.3");
end