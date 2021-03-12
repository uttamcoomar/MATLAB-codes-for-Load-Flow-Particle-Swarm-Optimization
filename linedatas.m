% Returns Line datas of the system...
function linedt = linedatas(num)
% | From | To | R | X | B/2 | X'mer |
% | Bus | Bus | pu | pu | pu | TAP (a) |
linedat5 = [ 1 2 .1 0.4 0 1
             1 4 0.15 0.6 0 1
             1 5 0.05 0.2 0 1
             2 3 0.05 0.2 0 1
             2 4 0.10 0.4 0 1
             3 5 0.05 0.2 0 1];
switch num
   case 5
  linedt = linedat5;
   case 30
  linedt = linedat30;
   case 57
  linedt = linedat57;
end
