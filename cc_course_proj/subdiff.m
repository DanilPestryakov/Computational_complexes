## -*- texinfo -*-
## @documentencoding UTF-8
## @defunx subdiff (@var{A}, @var{b})
## @defunx subdiff (@var{A}, @var{b}, @var{options})
##
## Returns formal solution of the interval linear system A &middot; x = b.
## in Kaucher arithmetics
##
## @param @var{A} the interval matrix A of the system
## @param @var{b} the right-hand side interval vector
## @param @var{options} parameters of subdifferential Newton's method 
## Default options:
##   options.max_iterations_num = 50
##   options.eps = 1e-6
##   options.tau = 1
##   options.norm_min_val = 1e-6
##
## @example
## @group
## Ainf = [ kinterval( 2,4) kinterval(-2,1);
##          kinterval(-1,2) kinterval( 2,4)];
## binf = [ kinterval(-2,2) kinterval(-2,2)];
## x = subdiff(A, b)
## @result{} x =
##  [-0.333333,0.333333]
##  [-0.333333,0.333333]
## Ax = A * x
## @result{} Ax =
##  [-2,2]
##  [-2,2]
## 
## opts.stepwise = 1; # solve the same ILS step by step
## [x, opts] = sdnewton1(A, b, opts);  # First iteration
## fprintf("=============\niteration: %d\n", opts.status.iteration);
## x 
##
## while ~opts.finish # Second and all the rest iterations
##    [x, opts] = sdnewton1(opts);
##    fprintf("=============\niteration: %d\n", opts.status.iteration);
##    x 
## endwhile
## @end group
## @end example
## @seealso{@@kinterval/resize, @@kinterval/cat, @@kinterval/postpad, @@kinterval/prepad}
## @end defunx

function [result, opt, iter, norm_diff] = subdiff(varargin)
global A b iteration n q r xx

    if nargin > 1 
        A = varargin{1};
        b = varargin{2}; 

        result = [];
        %if ~checkDimensions(A,b)
        %  return;
        %endif

        opt = setOptions(varargin);

        n = rows(A);

        F = superMatrix(A);
        xx = [b(:,:,1); b(:, :, 2)];
        xx = F \ xx;

        r = realmax;
        q = 1;

        iteration = 0;
    else
        opt = varargin{1};
        if (~isfield(opt,'stepwise') || ~isfield(opt,'finish') || ~isfield(opt,'finish'))
            error("Wrong usage (single argument is allowed in stepwise mode only)");
            print_usage();
            return;
        endif
        if opt.finish
            return;
        endif
        restoreStatus(opt.status);
    endif

    while r / q > opt.eps
        iteration = iteration + 1;

        r = 0;
        x = xx;
        F = zeros(2*n,2*n);

        for i = 1:n
          s = [0 0];
          for j = 1:n
            g = [A(i, j, 1) A(i, j, 2)];
            h = [x(j) x(j + n)];
            t = times(g, h); 
            s = s + t;
            F = calcSubgrad(F, i, j, g, h);
          endfor
          b_int = [b(i, 1, 1) b(i, 1, 2)];
          t = innerminus(s, b_int);

          xx(i) = t(1, 1);
          xx(i + n) = t(1, 2);

          r = r + max(abs(t(1, 1)), abs(t(1, 2)));
        endfor
        #disp(' ')
        xx =  F \ xx;
        xx = x -  xx .* opt.tau;
        q = norm(xx,1); # L1 norm
        if q <= opt.norm_min_val
            q = 1;
        endif

        if (isfield(opt,'stepwise') && opt.stepwise)
            opt.finish = ~(iteration < opt.max_iterations_num && r / q > opt.eps);
            opt.status = storeStatus();
            result = kinterval( xx(1:n), xx(n+1:end) );
            return;
        endif
    if iteration >= opt.max_iterations_num
      break;
    endif
    endwhile

    result = [xx(1:n); xx(n+1:end)];
    iter = iteration;
    norm_diff = r / q;
endfunction

function result = checkDimensions(A,b)
  if ~issquare(A)
      error ("matrix is not square");
      result = false;
  endif
  if columns(A) ~= rows(b)
      error ("mismatch of matrix and vector dimensions");
      result = false;
  endif
  result = true;
endfunction

function opt = setOptions(args)
    opt = struct("max_iterations_num",50, "eps",1e-6, "tau",1, "norm_min_val",1e-6, "mode","standard", "stepwise",false);

    if numel(args) > 2
        options = args{3};
        fn = fieldnames(opt);
        for k = 1:numel(fn)
            if isfield(options,fn{k})
                val = getfield(options, fn{k});
                opt = setfield (opt, fn{k}, val);
            endif
        endfor
    endif
endfunction

function result = superMatrix(A)
    midA = mid(A);
    P = (midA >= 0);
    midAp = zeros(size(midA));
    midAm = midAp;
    midAp(P) = midA(P);
    midAm(~P) = midA(~P);
    result = [midAp midAm; midAm midAp];
endfunction

function result = mid(A)
    result = zeros(rows(A), columns(A));
    for i = 1:rows(A)
      for j = 1:columns(A)
        result(i, j) = (A(i, j, 1) + A(i, j, 2)) / 2;
      endfor
    endfor
endfunction

function status = storeStatus()
    vars = who("global");

    eval(["global" sprintf(' %s',vars{:})]);

    vals = cell(size(vars));
    for k = 1:numel(vars)
        vals{k} = eval(vars{k});
    endfor

    status = cell2struct(vals, vars);
endfunction

function restoreStatus(status)
    vars = fieldnames(status);
    vals = struct2cell(status);

    eval(["global" sprintf(' %s',vars{:})]);

    for k = 1:numel(vars)
        assignin ("caller", vars{k}, vals{k}); 
    endfor
endfunction


# Set pair of F elements using A(i,j), [x(j),x(j+n)]
function F = calcSubgrad(F, i, j, a, b)
    n = rows(F) / 2;

    ai = a(1, 1);
    as = a(1, 2);
    bi = b(1, 1);
    bs = b(1, 2);
    
    if sign(ai) * sign(as) > 0 
         if sign(ai) > 0,  k = 0; else k = 2; endif
    else
         if ai < as,       k = 1; else k = 3; endif
    endif

    if sign(bi) * sign(bs) > 0
        if sign(bi) > 0,  m = 1; else m = 3; endif 
    else
         if bi <= bs,     m = 2; else m = 4; endif
    endif
    #disp(4 * k + m)
    switch (4 * k + m)
        case 1
            F(i,     j)     = ai;
            F(i + n, j + n) = as;
        case 2
            F(i,     j)     = as;
            F(i + n, j + n) = as;
        case 3
            F(i,     j)     = as;
            F(i + n, j + n) = ai;
        case 4
            F(i,     j)     = ai;
            F(i + n, j + n) = ai;
        case 5
            F(i,     j + n) = ai;
            F(i + n, j + n) = as;
        case 6
            if  ai * bs < as * bi
                F(i, j + n) = ai;
            else
                F(i, j)     = as;
            endif
            if ai * bi > as * bs
                F(i + n, j)     = ai;
            else
                F(i + n, j + n) = as;
            endif
        case 7
            F(i,     j) = as;
            F(i + n, j) = ai;
        case 9
            F(i, j + n) = ai;
            F(i + n, j) = as;
        case 10
            F(i, j + n) = ai;
            F(i + n, j) = ai;
        case 11
            F(i, j + n) = as;
            F(i + n, j) = ai;
        case 12
            F(i, j + n) = as;
            F(i + n, j) = as;
        case 13
            F(i,     j) = ai;
            F(i + n, j) = as;
        case 15
            F(i,     j + n) = as;
            F(i + n, j + n) = ai;
        case 16
            if ai * bi  > as * bs
                F(i, j)     = ai;
            else
                F(i, j + n) = -as;
            endif
            if ai * bs < as * bi
                F(i + n, j + n) = ai;
            else
                F(i + n, j)     = as;
            endif
    endswitch
endfunction

function z = times(x, y)
  xP = ispositive(x);
  yP = ispositive(y);
  xN = isnegative(x);
  yN = isnegative(y);
  xZ = iszeroin(x);
  yZ = iszeroin(y);
  if  xP 
    if yP
      z = [x(1, 1) * y(1, 1) x(1, 2) * y(1, 2)];
    elseif yZ 
      z = [x(1, 2) * y(1, 1) x(1, 2) * y(1, 2)];
    elseif yN
      z = [x(1, 2) * y(1, 1) x(1, 1) * y(1, 2)];
    else
      z = [x(1, 1) * y(1, 1) x(1, 1) * y(1, 2)];
  endif
  elseif xZ
     if yP
      z = [x(1, 1) * y(1, 2) x(1, 2) * y(1, 2)];
     elseif yZ 
       z = [min(x(1, 1) * y(1, 2), x(1, 2) * y(1, 1)) max(x(1, 1) * y(1, 1), x(1, 2) * y(1, 2))];
     elseif yN
       z = [x(1, 2) * y(1, 1), x(1, 1) * y(1, 1)];
     else
       z = [0 0];
  endif
  elseif xN
    if yP
      z = [x(1, 1) * y(1, 2) x(1, 2) * y(1, 1)];
    elseif yZ 
      z = [x(1, 1) * y(1, 2) x(1, 1) * y(1, 1)];
    elseif yN
      z = [x(1, 2) * y(1, 2) x(1, 1) * y(1, 1)];
    else
      z = [x(1, 2) * y(1, 2) x(1, 2) * y(1, 1)];
    endif
  else #if adualZ
    if yP
      z = [x(1, 1) * y(1, 1) x(1, 2) * y(1, 1)];
    elseif yZ 
      z = [0 0];
    elseif yN
      z = [x(1, 2) * y(1, 2) x(1, 1) * y(1, 2)];
    else
      z = [max(x(1, 1) * y(1, 1), x(1, 2) * y(1, 2)) min(x(1, 1) * y(1, 2), x(1, 2) * y(1, 1))];
    endif
  endif
endfunction
%
%    function z = inv(x)
%      if iszeroin(pro(x))
%        error ('wrong division, pro(denominator) contains zero');
%      else
%        z = kinterval(1/x.kinf, 1/x.ksup);	      
%      endif
%    endfunction
%
%    function z = rdivide(x,y)
%      if (not (isa (x, "kinterval")))
%        x = kinterval (x);
%      endif
%      if (not (isa (y, "kinterval")))
%        y = kinterval (y);
%      endif
%      if iszeroin(pro(y))
%        error ('wrong division, pro(denominator) contains zero');
%      endif
%      inv_y = kinterval(1/y.ksup, 1/y.kinf);
%      z = x .* inv_y;
%    endfunction
%
%    function z = mrdivide(x,y)
%      z = rdivide(x,y);
%    endfunction
%    
%    function z = opp(x)
%	    z = kinterval(-x.kinf, -x.ksup);
%    endfunction
%
function z = dual(x)
  z = [x(1, 2) x(1, 1)];
endfunction
%
%    function z = pro(x)
%      if x.kinf <= x.ksup
%          z = x;
%      else
%          z = kinterval(x.ksup, x.kinf);
%      endif
%    endfunction
%
function t = ispositive(x)
  t = x(1, 1) >= 0 &&  x(1, 2) >= 0;
endfunction
    
function t = iszeroin(x)
  t = x(1, 1) <= 0 &&  0 <= x(1, 2);
endfunction

function t = isnegative(x)
  t = x(1, 1) <= 0 &&  x(1, 2) <= 0;
endfunction

function t = isinzero(x)
  t = iszeroin(dual(x));
endfunction

function z = innerminus (x, y)
  zinf = x(1, 1) - y(1, 1);
  zsup = x(1, 2) - y(1, 2);

  z = [zinf zsup];
endfunction
