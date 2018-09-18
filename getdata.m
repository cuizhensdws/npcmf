function [Interaction,kCompound,kDisease,Did,Tid]=getdata(path,dataset)
%getdata loads the the interaction information of a particular dataset
%along with the pairwise similarities between the involved miRNA and
%Disease
%
% [Interaction,kCompound,kDisease,Did,Tid] = getdata(path,dataset)
%
% INPUT
%  path:        the path where all the files to be loaded are






    %------------------%
    % Adjacency Matrix %
    %------------------%

    % Extract the adjacency matrix...
    % The following command merges adjacency matrices that are 
    % extracted from multiple files which are specified using the 
    % pattern: [ path dataset '_admat_dgc.txt']
    newData1 = importdata ([ path dataset '_admat_dgc.txt']);
    Interaction = newData1.data;
    % Now, 'Interaction' has the matrix WITHOUT the column headers and
    % the row labels.


    Did = newData1.textdata(1,:);
    Did(1)=[];  % remove the first element (which is 'empty')


    Tid=newData1.textdata(:,1);
    Tid(1)=[];  % remove the first element (which is 'empty')



    newData1 = importdata ([ path dataset '_simmat_dc.txt']);
    kCompound = newData1.data;



    % Extract the disease similarity matrices and merge them into one.
    newData1 = importdata ([ path dataset '_simmat_dg.txt']);
    kDisease = newData1.data;



    clear newData1

    kCompound = (kCompound + kCompound')/2;

    epsilon = .1;
    while sum(eig(kCompound) >= 0) < size(kCompound,1) || isreal(eig(kCompound))==0
        kCompound = kCompound + epsilon*eye(size(kCompound,1));
    end
    while sum(eig(kDisease) >= 0) < size(kDisease,1) || isreal(eig(kDisease))==0
        kDisease = kDisease + epsilon*eye(size(kDisease,1));
    end
end