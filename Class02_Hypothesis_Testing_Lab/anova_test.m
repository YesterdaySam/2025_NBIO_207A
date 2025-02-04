% illustration of anovan function


% three groups, three separate times per group
g1{1} = normrnd(5 * ones(1,10), 3 * ones(1,10));
g1{2} = normrnd(5 * ones(1,10), 3 * ones(1,10));
g1{3} = normrnd(5 * ones(1,10), 3 * ones(1,10));

g2{1} = normrnd(5 * ones(1,10), 3 * ones(1,10));
g2{2} = normrnd(5 * ones(1,10), 3 * ones(1,10));
g2{3} = normrnd(5 * ones(1,10), 3 * ones(1,10));

g3{1} = normrnd(5 * ones(1,10), 3 * ones(1,10));
g3{2} = normrnd(5 * ones(1,10), 3 * ones(1,10));
g3{3} = normrnd(5 * ones(1,10), 3 * ones(1,10));

y = [];
group = [];
time = [];
for i = 1:3
    y = [y g1{i}];
    group = [group ones(size(g1{i}))];
    time = [time i*ones(size(g1{i}))];
end
for i = 1:3
    y = [y g2{i}];
    group = [group 2*ones(size(g2{i}))];
    time = [time i*ones(size(g2{i}))];
end
for i = 1:3
    y = [y g3{i}];
    group = [group 3*ones(size(g3{i}))];
    time = [time i*ones(size(g3{i}))];
end

[p atab stats] = anovan(y, {group time}, 3, 3, strvcat('group', 'time'));
comp = multcompare(stats, .05, 'off',[] ,[] ,[1 2])

