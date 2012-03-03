function p = ldl_search_direction(H,g,delta)

[L,D,P] = ldl(H);
% compute D_tilde
[row,column,s] = find(D);
off_diags = find(row-column<0);

just_diag_entries = ones(size(D,1),1);

if length(off_diags)>0
for off_diag_entry = off_diags
  % zero out the entries associated with the off diagon
  just_diag_entries(row(off_diag_entry)) = 0;
  just_diag_entries(column(off_diag_entry)) = 0;

  % entry corresponds to upper right hand corner
  r = row(off_diag_entry);
  
  c = column(off_diag_entry)-1;
  
  D_mini_block = D(r:r+1,c:c+1);
  
  
  [V,D_tilde] = eigs(D_mini_block);
  
  D_tilde = D_tilde + max(diag([delta ; delta]) -D_tilde,0);

  D(r:r+1,c:c+1) = V * D_tilde_thresh *V';
end

  % threshold all the just diagonal guys to make sure they are large enough
end
  D_diag = diag(D);
  D_posdef = D + diag( just_diag_entries .* max(delta *ones(size(D,1),1) - D_diag,0));
  %display(num2str(min(diag(D_posdef))));



DLtp = L\g;
Ltp = D_posdef \ DLtp;
p = L' \ Ltp;


%LSopts.SYM = true;
%LSopts.POSDEF = true;

%Lopts.LT=true;
%DLtp = linsolve(L,g,Lopts);
%Lopts.LT=false;
%Lopts.POSDEF = true;
%Lopts.SYM = true;
%Ltp = linsolve (D_posdef,DLtp,Lopts);
%Lopts.UT = true;
%Lopts.POSDEF = false;
%Lopts.SYM = false;
%p = linsolve (L',Ltp,Lopts);

