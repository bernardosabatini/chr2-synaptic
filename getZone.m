function zone=getZone(newCell)

% use this for the CTB defined zone	
zone=newCell.Injection(1);
	
% use this for the anatomically defined zone
if newCell.ML<=0
	zone='M';
elseif newCell.ML<=400
	zone='C';
else 
	zone='L';
end		