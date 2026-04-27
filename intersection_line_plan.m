function [M,TST]=intersection_line_plan(Ap,Bp,Cp,Al,Bl,eps)
    % Ap,Bp,Cp : three points defining the plan
    % Al,Bl : two points defining the line
    % eps : round-off error tolerance
    % M : intersection point
    % TST: Test variable: true if M is inside (edges included), false if not.
    
    abc = cross(Bp-Ap, Cp-Ap);
        N=abc/norm(abc); % Normal to the plan
    d = -sum(Ap.*N); 
    
    U=Bl-Al; U=U/norm(U); % Line vector
    
    div=sum(N.*U);
    % Colinearity test
    if abs(div)<10^(-8) % To avoid decimal error that could happen
    	M=[]; TST=false;
    else
        % Plan equation: sum(N.*M)+d=0
        % Line equation: M-Al=t*U
        t=-(sum(N.*Al)+d)/div;
        M=t*U+Al;
        
        % Test if M is inside the triangle
        Tri_tst=inside_triangle(Ap,Bp,Cp,M,eps);
        % Test if M is between Al and Bl
        Lin_tst=(all(M >=(min([Al;Bl])-eps)) && all(M <=(max([Al;Bl])+eps)));
        
        TST = Tri_tst && Lin_tst;
    end
end
