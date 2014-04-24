function submit = Midway(core,mode)
submit =['module load matlab; '];
for i = 0:mode
    submit = [submit 'sbatch ' core '_' num2str(i) '.sbatch; '];
    fid = fopen([core '_' num2str(i) '.sbatch'],'W');
    fprintf(fid,'#!/bin/bash\n');
    fprintf(fid,'#SBATCH --output=%s\n',[core '_' num2str(i) '.txt']);
    fprintf(fid,'#SBATCH --job-name=%s\n',[core '_' num2str(i)]);
    fprintf(fid,'#SBATCH --partition=bigmem\n');
    % fprintf(fid,'cd $Home/rcchelp/software/matlab.rcc-docs');
    fprintf(fid,'module load matlab\n');
    fprintf(fid,'matlab -nodisplay  -r "%s(%d)"\n', core,i);
    fclose(fid);
end
disp(submit);
return