allFields=fieldnames(apResults.all);



for fc=1:length(allFields)
	fns=allFields{fc};
	
	disp(fns)
	for cc={'all', 'M', 'C', 'L'}
		[mean(apResults.(cc{1}).(fns)) std(apResults.(cc{1}).(fns))]
	end
	disp('')
end