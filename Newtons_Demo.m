% Newton's Method - Find intersections of two circles
% 
% (x1-x0)^2 + (y1-y0)^2  = r1^2
% (x2-x0)^2 + (y2-y0)^2  = r2^2
%
% Author: David Ohm
% -----------------------------------------------------

% Starting coordinates
x = 7.2;
y = 10.1;

% Start iteration
for i = 1:5
   
    X = [x y];
    % Define two overlapping circles
    r1 = 4;             % radius of circles
    r2 = 4;
    x1 = 1; y1 = 4;     % center of circles
    x2 = 3; y2 = 2;

    % Compute Jacobian Matrix
    A11 = 2*(x-x1);
    A12 = 2*(y-y1);
    A21 = 2*(x-x2);
    A22 = 2*(y-y2);
    A = [A11 A12;A21 A22];

    % Compute functions
    F1 = r1.^2 - (x-x1).^2 - (y-y1).^2;
    F2 = r2.^2 - (x-x2).^2 - (y-y2).^2;

    % Compute Xnew
    B = [F1 F2];
    Ainv = inv(A);
    DeltaX = Ainv*B';
    Xnew = X+DeltaX'
    
    x = Xnew(1);
    y = Xnew(2);
    
end 


