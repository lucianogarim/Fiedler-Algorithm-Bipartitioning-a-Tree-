Tree_Generator_main

[autv,autva]=eigs(laplacian,size(laplacian,1));
U2=autv(:,size(laplacian,1)-1);

%Fiedler's algorithm (RatioCut formulation)
for j=1:size(laplacian,1)
    t(j)=U2(j,1);                       %set the threshold according to the vector components
    x = (U2 > t(j));                    %taking the elements higher than the threshold
    x = (2 * x) - 1;                    %adjusting the vector x to values ??of 1 and -1 (indicators)
    s=sum(x>-1);                        %defining the sizes of the sets
    scomp=sum(x < 1);
    n=s+scomp;
    for i=1:n
        if x(i)>0
            x(i)= sqrt((scomp)/(s*n));  %building the characteristic vectors
        else
            x(i)=-sqrt((s)/(scomp*n));
        end
    end
    rrcut(j) = (x' * laplacian * x);    %calculating the ratiocut %%% The generated partition is in vector x
end
rrcut=transpose(rrcut);

[M,I] = min(rrcut(rrcut>0));
t(I)=U2(I,1);
x = (U2 > t(I));
x = (2 * x) - 1;
s=sum(x>-1);
scomp=sum(x < 1);
n=s+scomp;
for i=1:n
    if x(i)>0
        x(i)= sqrt((scomp)/(s*n));
    else
        x(i)=-sqrt((s)/(scomp*n));
    end
end
rrcut(I) = (x' * laplacian * x);

G=graph(adj);
p = plot(G);                        %Plot the result
b=zeros(n,1);
for i=1:n
    if x(i)>=0
        b(i)=i;
    end
end
b = nonzeros(b');
highlight(p,b,'NodeColor','r')