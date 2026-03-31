function [var_estimate] = ESPER_wrapper_for_python(DesiredVariables,Coordinates_all,Measurements_all,MeasIDVec_ESPER,Equations , Dates, ESPER_type, VerboseTF, pHCalcTF)

Dates = Dates'; % cludge, but seems to work for now 

% Testing:
% DesiredVariables = 1;
% Coordinates_all = [-145, 50, 0];
% Measurements_all = [34, 10, 300];
% MeasIDVec_ESPER = [1, 2, 6];
% Equations = 7;
% Dates = 2020.1;
% VerboseTF = 0;
% pHCalcTF = 0;

% disp(DesiredVariables)
% disp(MeasIDVec_ESPER)
% disp(Equations)
% disp(VerboseTF)
% 
% disp(size(Coordinates_all))
% disp(size(Measurements_all))
% disp(size(Dates))
if strcmp(ESPER_type,'MX')
    [Estimates,Uncertainties]=ESPER_Mixed(DesiredVariables,Coordinates_all,Measurements_all,MeasIDVec_ESPER,'Equations', Equations, 'EstDates', Dates, 'VerboseTF', VerboseTF, 'pHCalcTF', pHCalcTF);
elseif strcmp(ESPER_type,'NN')
    [Estimates,Uncertainties]=ESPER_NN(DesiredVariables,Coordinates_all,Measurements_all,MeasIDVec_ESPER,'Equations', Equations, 'EstDates', Dates, 'VerboseTF', VerboseTF, 'pHCalcTF', pHCalcTF);
elseif strcmp(ESPER_type,'LIR')
    [Estimates,Uncertainties]=ESPER_LIR(DesiredVariables,Coordinates_all,Measurements_all,MeasIDVec_ESPER,'Equations', Equations, 'EstDates', Dates, 'VerboseTF', VerboseTF, 'pHCalcTF', pHCalcTF);
end
% disp(Estimates)

var = fieldnames(Estimates);
var_estimate = Estimates.(var{1});
var_uncertainty = Uncertainties.(var{1});

end