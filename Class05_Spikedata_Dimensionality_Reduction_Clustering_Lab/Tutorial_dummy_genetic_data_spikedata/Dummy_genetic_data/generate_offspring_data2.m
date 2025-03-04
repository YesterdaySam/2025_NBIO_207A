% offspring_data.m
rng(1)

Nsites = 4000;
Ntypes = 4;
N2 = 4;
N3 = randi(5,N2,1);
N3tot = sum(N3);
N4 = randi(5,N3tot,1);
N4tot = sum(N4);

mutate_rate = 1/500;

g1 = randi(Ntypes,Nsites,2);

g2 = zeros(Nsites,N2);
type2 = randi(2,Nsites,N2);

for i = 1:N2
    for j = 1:Nsites
        g2(j,i) = g1(j,type2(j,i));
    end
end

g2 = g2 + (rand(Nsites,N2) < mutate_rate).*(randi(Ntypes,Nsites,N2)-g2);

g3 = zeros(Nsites,N3tot);
type3 = randi(2,Nsites,N3tot) -1;

k = 0;
for i = 1:N2
    for j = 1:N3(i)
        k = k+1;
        indices = find(type3(:,k));
        g3(indices,k) = g2(indices,i);
        indices = find(~type3(:,k));
        g3(indices,k) = randi(4,length(indices),1);
    end
end

g3 = g3 + (rand(Nsites,N3tot) < mutate_rate).*(randi(Ntypes,Nsites,N3tot)-g3);

g4 = zeros(Nsites,N4tot);
type4 = randi(2,Nsites,N4tot)-1;

k = 0;
for i = 1:N3tot
    for j = 1:N4(i)
        k = k+1;
        indices = find(type4(:,k));
        g4(indices,k) = g3(indices,i);
        indices = find(~type4(:,k));
        g4(indices,k) = randi(4,length(indices),1);
    end
end


g4 = g4 + (rand(Nsites,N4tot) < mutate_rate).*(randi(Ntypes,Nsites,N4tot)-g4);

g_tot = [g1 g2 g3 g4];
N_tot = 2 + N2 + N3tot + N4tot;

g4_new = zeros(N4tot,Nsites*Ntypes);
for i = 1:Nsites
    start = (i-1)*Ntypes;
    for j = 1:Ntypes
        indices = find(g4(i,:) == j );
        g4_new(indices,start+j) = 1;
    end
end

Z = linkage(g4_new,'average','squaredeuclidean');
figure(5)
dendrogram(Z,N_tot)


        

