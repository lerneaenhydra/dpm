function rads = rpm2rads(rpm)
%RPM2RADS rads = rpm2rads(rpm) scales an input representing an angular
%speed in revolutions per minute to radians/second

rads = rpm * 2 * pi / 60;

end

