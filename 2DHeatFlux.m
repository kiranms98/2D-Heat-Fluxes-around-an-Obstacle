Coorde=importdata('FineNodeCoord.txt');
Conne=importdata('FineNodeConnRound.txt');
Conn=Conne(:,2:4);
Coord=Coorde(:,2:3);
x=Coord(:,1);
y=Coord(:,2);
k=5;
K=AssembleStiffness(Conn, Coord,k);
t=[-1/sqrt(3);1/sqrt(3)];
F_elem(Coord,Conn(1,:),-1,t);
Connectivity_element=Conn(1,:);
F_tot=Ftot(Coord,Conn,-1,t);
[K_reduc,F_reduc]=Biff(K,F_tot,Conn,Coord);
T= SolveT(K_reduc,F_reduc);
T_tot= ComputeTotalT(T,Conn,Coord);
[q]= Computeq(Conn,Coord,T_tot,k);
%% Definition of function for number of equation
function numEq= ComputeNumEq(Connectivity)
    n_elem  = size(Connectivity,1); % The number of elements
    n_nodes_per_elem = size(Connectivity,2);
    numEq = zeros(n_elem, n_nodes_per_elem);
    for e=1:n_elem
        for i=1:n_nodes_per_elem
            numEq(e, i) = Connectivity(e, i);
            
        end
    end
end
%% Definition of function to compute gradients
function C = ComputeC()
    N1s=-1/2;
    N1t=-1/2;
    N2s=1/2;
    N2t=0;
    N3s=0;
    N3t=1/2;
    C = [
        [N1s,N2s,N3s];
        [N1t,N2t, N3t]];
        
end

%% Definition of function to compute B and J
function [B,J]= ComputeBandJ(~, nodes)
    C = ComputeC();
    J = C*nodes;
    Mat_J = inv(J);
    B = Mat_J*C; 
end
%% Definition of local stiffness matrix
function [K_local]= ComputeLocalMatrix(connectivity_element, coord,k)
    
    quads = [[1/sqrt(3), -1/sqrt(3)]];
    weights = [[1,1]];
    idx=reshape(connectivity_element,[3,1]);
    nodes = coord(idx,:); 
    K_local = zeros(3, 3); 
    
    for i=1:size(quads)
        for j=1:size(weights)
        [B, J] = ComputeBandJ(quads(i), nodes);  
        detJ = det(J);          
        K_local = K_local + k*weights(j) * B'*B * detJ;
        end
    end
end
%% Compute global stiffness matrix 
function [K_global]= AssembleStiffness(Connectivity, Coord,k)

    n_elem  = size(Connectivity,1);
    n_nodes = size(Coord,1);
    numEq = ComputeNumEq(Connectivity);

    K_global = zeros(n_nodes, n_nodes);

    for e=1:n_elem 
        Degree_of_freedom = numEq(e, :); %Equation associated to the degree of freedom
        Connectivity_element = Connectivity(e, :); % Connectivity of the element: which nodes are we studying?
        K_local = ComputeLocalMatrix(Connectivity_element, Coord,k);  %Compute local matrix
        
        for i=1:size(Degree_of_freedom,2)
            for j=1:size(Degree_of_freedom,2) 
                gi=Degree_of_freedom(1,i);
                gj=Degree_of_freedom(1,j);
                K_global(gi,gj)=K_global(gi, gj)+ K_local(i, j);
            end
        end
    end
end

 %%   Compute source matrix
 function [N_fct]= N(s,t)
    N_fct= [
    [ - s/2 - t/2];
    [s/2+1/2];
    [t/2+1/2]];
 end
 function [F_elem]= F_elem(Coord,Connectivity_element,q,t)
    
    idx=reshape(Connectivity_element,[3,1]);
    nodes = Coord(idx,:);
    y3=nodes(3,2);
    y1=nodes(1,2);
    J=abs(y3-y1)/2;
    F_elem=zeros(3,1);
    for i=1:size(t)
        [N_fct]=N(-1,t(i));
        F_elem= F_elem+-q*[[J*N_fct(1)];[0];[J*N_fct(3)]];
    end
 end
 
 
 function [F_tot]=Ftot(Coord,Connectivity,q,t)
         n_elem  = size(Connectivity,1);
         n_nodes = size(Coord,1);
         numEq = ComputeNumEq(Connectivity);
         xF=Coord(:,1)
         yF=Coord(:,2)
         F_tot=zeros(n_nodes,1);
        
        for e=1:n_elem
            
            Connectivity_element=Connectivity(e, :); %Get the connectivity of the element we are considering
            F_element=F_elem(Coord,Connectivity_element,q,t); %Compute the element flux of the element we are considering
            Degree_of_freedom = numEq(e, :);
            for i=1:size(Degree_of_freedom,2)
                if xF(Degree_of_freedom(i))==-10
                        gi=Degree_of_freedom(i);
                        F_tot(gi)=F_tot(gi)+F_element(i);
                end
            end
        end
 end   
            
 function [K_reduc,F_reduc]=Biff(K,Ftot,Conn,Coord) %We delete all the rows associated with a T=0 (the blocked ones)
    n_elem  = size(Conn,1); % The number of elements
    n_nodes=size(Coord,1)
    xF=Coord(:,1);
    K_reduc=K; %Initialise the reduced K with K
    F_reduc=Ftot %Initialize the reduced F with Ftot
    %for e=1:n_elem
        
            %Degree_of_freedom = Conn(e, :);
            %for i=1:size(Degree_of_freedom,2)
                %gi=Degree_of_freedom(i)
       for i=1:n_nodes
          if xF(i)==10 %If we are at the end of the bar--> T=0 
            K_reduc(i,:)=[]; %DELETE row associated with blocked node
            K_reduc(:,i)=[]; %Delete column associated with blocked node
            F_reduc(i)=[]; %Delete row associated with blocked node, we can do that bc we know we don't have any source and the only
                                         %applied fluxes are on the
                                         %leftmost part of the beam
                end
            end
    end
    
   
                

 function T= SolveT(K,F)
    T=inv(K)*F;
 end
 function T_tot= ComputeTotalT(T,Conn,Coord)
    
    n_nodes = size(Coord,1);
    T_tot=zeros(n_nodes,1);
    c=1;
    xc=Coord(:,1);
    for i=1:n_nodes
        if(xc(i)==10)
            T_tot(i)=0;
        else
           
           T_tot(i)=T(c);
           c=c+1;
           
        end
    end
 end
 function [q]= Computeq(Conn,Coord,T,k)
    n_elem  = size(Conn,1);
    n_nodes = size(Coord,1);
    numEq = ComputeNumEq(Conn);
    Connectivity_element = Conn(1, :)
    idx=reshape(Connectivity_element,[3,1]);
    q=zeros(2*n_elem,1);
    nodes = Coord(idx,:)
    [B,J]= ComputeBandJ(0, nodes);
    q(1:2)=-k*B*T(idx); %First two rows of the flux: flux of the first element
    
    for e=2:2:n_elem
        Connectivity_element = Conn(e, :); %Look for connectivity of the element we are considering to know which nodes we are talking about
        idx=reshape(Connectivity_element,[3,1]);
        nodes = Coord(idx,:) %Get the coordinates of these nodes
        [B,J]= ComputeBandJ(0, nodes); %Compute the B and the Jacobian of this element
        q(e:e+1)=-k*B*T(idx); %Compute the flux of this element
        %Repeat but for another element and each time add the two flux
        %values obtained at the righ place of Q ( position=Node_number)
    end
 end
