function [AmExp, OrcaExp] = autoExposure(objCCM, objRef, AmExp, OrcaExp)

% define constants for the upper and lower bounds of exposure times
min_orca_exp = 0.01;
max_orca_exp = 3;
min_am_exp = 0.01;
max_am_exp = 4000;

% define target exposures
target_ccm = 40000;
target_ref = 200;

% adjust OrcaExp based on objCCM
max_ccm = max(objCCM, [], 'all');
adjustment_orca = (max_ccm - target_ccm) * 0.0001;  % calculate proportional adjustment, adjust this constant as needed
new_orca_exp = OrcaExp - adjustment_orca;  % calculate new proposed exposure time

% ensure new exposure time is within the specified bounds
OrcaExp = min(max(new_orca_exp, min_orca_exp), max_orca_exp);

% adjust AmExp based on objRef
max_ref = max(objRef, [], 'all');
adjustment_am = (max_ref - target_ref) * 0.5;  % calculate proportional adjustment, adjust this constant as needed
new_am_exp = AmExp - adjustment_am;  % calculate new proposed exposure time

% ensure new exposure time is within the specified bounds
AmExp = min(max(new_am_exp, min_am_exp), max_am_exp);
