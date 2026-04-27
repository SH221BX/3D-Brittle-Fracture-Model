function TST=inside_triangle(A,B,C,M,eps)
% Test if M is inside (edges included) the triangle formed by the three
% points A, B, C.
% M must already be in the plan ABC or it doesn't work
% eps : round-off error tolerance

% In 2D, it is easy to check this: the cross product AMB, BMC, CMA must be
% of identical signs.
% To extend this idea in 3D, I project the triangle on the three plans 
% XY, YZ and ZX and does the 2D test (I migh be wrong but it seems to work
% really well).

iL=[1 2; 2 3; 3 1]; 
TST= true;
for l=1:3    
    p=iL(l,:); a=A(p); b=B(p); c=C(p); m=M(p);
    ma=(a-m)/vecnorm(a-m); mb=(b-m)/vecnorm(b-m); mc=(c-m)/vecnorm(c-m); 
    ca=(a-c)/vecnorm(a-c); cb=(b-c)/vecnorm(b-c);

    % Test if a,b,c are not colinear
    v_tst=ca(1)*cb(2)-ca(2)*cb(1);
    if abs(v_tst) > eps % To avoid decimal error that could happen 
        v=[ ma(1)*mb(2)-ma(2)*mb(1)
            mb(1)*mc(2)-mb(2)*mc(1)
            mc(1)*ma(2)-mc(2)*ma(1) ];
        v(abs(v)<eps)=0;
        TST = TST && ( ...
                    (sign(v(1)*v(2))==1 || sign(v(1)*v(2))==0) ...
                    && ...
                    (sign(v(1)*v(3))==1 || sign(v(1)*v(3))==0) ...
                    && ...
                    (sign(v(2)*v(3))==1 || sign(v(2)*v(3))==0) ...
                    );
    else % If the projection forms a line
        TST = TST && ( ...    
              (all(m >=(min([a;b;c])-eps)) && all(m <=(max([a;b;c])+eps))) ...
                    );
    end
end

end