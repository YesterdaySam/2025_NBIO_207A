% Probability_Problem2
% -----------

% set the number of runs for the simulation
nruns = 10000;
% the probability that a fired is 0.1
pA = 0.1;
% the probability that b fired is 0.4
pB = 0.4;
% the probability of c firing given that a fired is 0.5
pCgivenA = 0.5;
% the probability of c firing given that b fired is 0.2
pCgivenB = 0.2;
% the probability of c firing given that a and b fired is 1
pCgivenAB = 1;
%create empty matricies to hold the zeros and ones corresponding to whether
%neurons A, B and C spike or not on each run
A = zeros(nruns,1);
B = zeros(nruns,1);
C = zeros(nruns,1);

% use a for loop to iterate through the runs for neurons a and b
% --------------------------------------------------------------
for run = 1:nruns
    % generate a random number bewteen 0 and 1
    randnum = rand(1,1);
    % set a to have fired if that number is < pa;
    if (randnum < pA)
        A(run) = 1;
    else
        A(run) = 0;
        % note that this is redundant, as a(run) is 0 to begin with;
    end
    % repeat that for neuron b
    randnum = rand(1,1);
    if (randnum < pB)
        B(run) = 1;
    else
        B(run) = 0;
    end
    % now that we have the firing of a and b, we can generate the simulations
    % for c. We'll do this with more compact code
    % c | a
    if (A(run) == 1)
        % compare the random number to pcgivena and assign the result
        C(run) = rand(1,1) < pCgivenA;
    end
    % if a did not cause c to fire, b still could
    if ((C(run) == 0) & (B(run) == 1))
        % compare the random number to pcgivenb and assign the result
        C(run) = rand(1,1) < pCgivenB;
    end
    % if neither a nor b caused c to fire on their own, but both were active, c
    % could fire
    if ((C(run) == 0) & (A(run) == 1) & (B(run) == 1))
        C(run) = rand(1,1) < pCgivenAB;
    end
end

% now we can answer the questions
% -----------------------------
% we'll go through all of the runs and count up the number where C fired and
% the number of those where A also fired, and similarly for B and the other
% combinations

% Initialize the variables and go through each simulation (1 to nrun)
% that was run above, and count how many times the conditions for each
% variable are met
% ---------------------------------------------------------
Cfired = 0; % c
ACfired = 0;% a & c
BCfired = 0; % b & c
ABCfired = 0; %a, b and c
AnBCfired = 0; % a, not b, and c
nABCfired = 0; % not a, b, and c
for run = 1:nruns
    % check to see if c fired on this run.
    % Note that if (a) is equivalent to if (a ~= 0)Prob
    if (C(run))
        Cfired = Cfired + 1;
        if (A(run))
            ACfired = ACfired + 1;
        end
        if (B(run))
            BCfired = BCfired + 1;
        end
        if (A(run) & B(run))
            ABCfired = ABCfired + 1;
        end
        if (A(run) & ~B(run))
            AnBCfired = AnBCfired + 1;
        end
        if (~A(run) & B(run))
            nABCfired = nABCfired + 1;
        end
    end
end
% p(A|C) is the proportion of times that a fired when c fired:
pAgivenC = ACfired / Cfired
0.4997
% p(B|C) is the proportion of times that a fired when c fired:
pBgivenC = BCfired / Cfired
0.7784
% p(AB|C) is the proportion of times that a and b fired when c fired:
pABgivenC = ABCfired / Cfired
0.278
% p(A & ~B|C) is the proportion of times that a but not b fired when c fired:
pAnBgivenC = AnBCfired / Cfired
0.2216
% p(~A & B|C) is the proportion of times that b but not a fired when c fired:
pnABgivenC = nABCfired / Cfired
0.5003
%



% ------------------------------------------------------------------------------


% That all works, but it's hideously inefficient. Here is a much more efficient
% version
nruns = 10000;
a = zeros(nruns,1);
b = zeros(nruns,1);
% a and not b and c
anbc = zeros(nruns,1);
% not a and b and c
nabc = zeros(nruns,1);
% a and b and c
abc = zeros(nruns,1);
f = rand(nruns,4);
% use the first two columns to figure out when a and b fired.
a = f(:,1) < pA;
b = f(:,2) < pB;
% to get anbc, we use column 3 to figure out when C fired as a result of A
% alone
anbc = a & ~b & (f(:,3) < pCgivenA);
% to get nanc, we use column 4 to figure out when C fired as a result of B
% alone
nabc = ~a & b & (f(:,4) < pCgivenB);
% to get a & b & c we only need to get a & b, as P(C | A & B) = 1
abc = a & b;
% the set of runs where C fired is the union of the anbc, nabd and abc arrays
c = anbc | nabc | abc;
% p(a|c)
pac = (sum(anbc) + sum(abc)) / sum(c)
% p(b|c)
pbc = (sum(nabc) + sum(abc)) / sum(c)
%p(ab|c)
pabc = sum(abc) / sum(c)
%p(a~b|c)
panbc = sum(anbc) / sum(c)
%p(~ab|c)
pnabc = sum(nabc) / sum(c)





% ------------------------------------------------------------------------------


% C. Error calculation. The most efficient way to do this is to do one run of
% 100000 iterations and take parts of it for each comparison.
% here is the answer from part A
actualpac = 0.493;
nruns = 100000;
a = zeros(nruns,1);
b = zeros(nruns,1);
% a and not b and c
anbc = zeros(nruns,1);
% not a and b and c
nabc = zeros(nruns,1);
% a and b and c
abc = zeros(nruns,1);
f = rand(nruns,4);
% use the first two columns to figure out when a and b fired.
a = f(:,1) < pA;
b = f(:,2) < pB;
anbc = a & ~b & (f(:,3) < pCgivenA);
nabc = ~a & b & (f(:,4) < pCgivenB);
abc = a & b;
c = anbc | nabc | abc;
nruns = [1 10 100 1000 10000 100000];
err = zeros(size(nruns));
for r = 1:length(nruns)
    err(r) = abs(actualpac - (sum(anbc(1:nruns(r))) + sum(abc(1:nruns(r)))) /...
        sum(c(1:nruns(r))));
end
% in our case, neuron C didn't fire on the first run, so err(1) is NaN.
% plot the errors that make sense
figure;
semilogx(nruns(2:end), err(2:end), 'LineWidth', 2)
h = xlabel('Number of runs')
set(h, 'FontSize', 20);
h = ylabel('Error');
set(h, 'FontSize', 20);
set(gca, 'YLim', [0 max(err(2:end))*1.1]);
set(gca, 'FontSize', 18);


% ------------------------------------------------------------------------------

 
% for x=1:1000; R(x)=rand(1,1); end
% figure; hold on; plot(R)
% figure; hold on; plot(sort(R))
% figure; hold on; plot(cumsum(sort(R)))




