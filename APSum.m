allFields=fieldnames(apResults.all);



for fc=1:length(allFields)
	fns=allFields{fc};
	
	disp(fns)
	c=1;
	ss=zeros(4,2);
	for cc={'all', 'M', 'C', 'L'}
		ss(c,:)=[mean(apResults.(cc{1}).(fns)) std(apResults.(cc{1}).(fns))];
		disp(ss(c,:))
		c=c+1;
	end
	disp('')
	
	figure
	hold on
	bar(1:4,ss(:,1))
	errorbar(1:4,ss(:,1),ss(:,2),'.')
	title(fns)

end

