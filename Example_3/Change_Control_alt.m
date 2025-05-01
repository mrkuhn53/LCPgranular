function Change_Control_Period = ...
         Change_Control_alt(istep, igoal, increment, max_increments, ...
                            Time, lap_Time, start_Time);
%
% this function determines whether the current control period has finished,
% so that the next control period should be entered
%
  Change_Control_Period = 0;
%
  if igoal(istep) == 70
    Change_Control_Period = increment >= max_increments;
  elseif igoal(istep) == 80
    Change_Control_Period = Time >= start_Time + lap_Time;
  end
%
  Change_Control_Period;
