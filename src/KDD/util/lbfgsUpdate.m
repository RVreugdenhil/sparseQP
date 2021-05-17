function [old_dirs,old_stps,Hdiag] = lbfgsUpdate(y,s,old_dirs,old_stps,Hdiag)
corrections=10;
debug=0;


ys = trace(y'*s);
if ys > 1e-10,
    numCorrections = size(old_dirs,3);%看有多少列
    if numCorrections < corrections
        % Full Update
        old_dirs(:,:,numCorrections+1) = s;
        old_stps(:,:,numCorrections+1) = y;
    else
        % Limited-Memory Update
        
        for i=1:corrections-1,
            old_dirs(:,:,i)=old_dirs(:,:,i+1);
            old_stps(:,:,i)=old_stps(:,:,i+1);
            
        end;
        old_dirs(:,:,corrections)=s;
        old_stps(:,:,corrections)=y;
    end
    
    % Update scale of initial Hessian approximation
    Hdiag = ys/trace(y'*y);
else
    if(debug)
        fprintf('Skipping Update\n');
    end
    
end
