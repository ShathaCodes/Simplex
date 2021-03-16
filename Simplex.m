function [k] = Simplex(C,A,v,b,isMax)

% C matrice contenant les coeff des var dans la fonction
% A matrice contenant les coeff des var dans les contraintes
% v matrice contenant le type de contrainte ( vi(i)=1 si <=, v(i)=-1 si >=)
% b matrice contenant les valeurs des contraintes
% isMax =1 si sonction de maximisation, isMax=0 sinon
%
% exemple :
% C=[1000 1200];
% A=[10 5; 2 3; 1 0; 0 1 ];
% v=[1;1;1;1];
% b=[200;60;34;14];
% isMax=1

%-------------------------------------------------------Definir variables

n=length(C);           % nbr variables
m=length(b);           % nbr equations
A=[A,diag(v,0),b];     % turn A into full matrice with ei and bj
C=[C,zeros(1,m)];      % add coeff ei
Z=zeros(1,n+m+1);      % initiallement des zeros
CZ=C-Z(1:n+m);
coeffb=zeros(m+2,1);   % coeff des var base
b=zeros(m+2,1);        % les variables de base (ex 1->X1 )

%-------------------------------------------------------affichage

M=[ [0 ; coeffb], [0 ; b], [ [C , 0] ; A ; Z ; [CZ,0] ] ];
disp(M);

%-------------------------------------------------------boucle jusqu'a fin

fin=1;
while ( fin == 1)
    
%-------------------------------------------------------calcul ve

	ve=1;
    for i=2:1:n+m
        if(isMax ==1 )
            if(CZ(i)>CZ(ve))
                ve=i;
            end;
        else
            if (CZ(i)<CZ(ve))
                ve=i;
            end
        end
    end
    if (ve <=n)
        fprintf('variable entrant ve= X%d\n',ve);
    else
        fprintf('variable entrant ve= e%d\n',ve-n);
    end
    
%-------------------------------------------------------calcul vs

    mn=A(1,n+m+1)/A(1,ve);
    vs=1;
    for i=2:1:m
        if( A(i,n+m+1)/abs(A(i,ve))<mn)
            vs=i;
            mn=A(i,n+m+1)/A(i,ve);
        end
    end
    if(coeffb(vs) ==0)
        fprintf('variable sortant vs= e%d\n',vs);
    else
        fprintf('variable sortant vs= X%d\n',b(vs));
    end
    %disp('vs=');
    %disp(vs);
    
%-------------------------------------------------------pivot

    for i=1:1:m
        for j=1:1:n+m+1
            if ( i~=vs && j~=ve && A(vs,ve)~=0)
                A(i,j) = A(i,j)-(A(i,ve)*A(vs,j)/A(vs,ve));
            end
        end
    end
%colonne of pivot : set to 0
    for i=1:1:m
        if(i~=vs)
            A(i,ve)=0;
        end
    end
%line of pivot : divise by itself so pivot=1
    pivot=A(vs,ve);
    for j=1:1:n+m+1
        A(vs,j) = A(vs,j)/pivot;
    end
    
%-------------------------------------------------------calcul coeffb et Z et CZ
    
    %set ve in place of vs
    coeffb(vs)=C(ve);
    b(vs)=ve;
    % calculate Z
    for j=1:1:n+m+1
        somme=0;
        for i=1:1:m
            somme = somme + coeffb(i)*A(i,j);
        end
        Z(j)=somme;
    end
    CZ=C-Z(1:n+m);

%-------------------------------------------------------affichage

    M=[ [0 ; coeffb], [0 ; b], [ [C , 0] ; A ; Z ; [CZ,0] ] ];
    disp(M);
    
%-------------------------------------------------------tester si fin

    fin=0;
    for i=1:1:n+m
        if((CZ(i) >0 && isMax==1) || (CZ(i) <0 && isMax==0) )
            fin=1;
        end
    end
end

%-------------------------------------------------------afficher Z et Xi

fprintf('Solution :\n');
k=Z(n+m+1);
for j=1:1:n
    vx=0;
    for i=1:1:m
        if (b(i) == j)
            vx=A(i,n+m+1);
        end
    end
    fprintf('X%d= %d\n',j,vx);
end
fprintf('Z= %d\n',k);
end
