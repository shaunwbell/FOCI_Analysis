#!/usr/bin/env python

import fileinput

print "<gpx version=\"1.1\" creator=\"Dave Strausz Python Program\">"

for line in fileinput.input() :
    line = line.rstrip().split("*")
    station = line[0]
    mooring = line[1]
    mooring_site = line[2]
    operation = line[3]
    depth = line[4]
    lat = line[5]
    lon = line[6]
    print "    <wpt lat=\"{0}\" lon=\"{1}\">" .format(lat, lon)
    if mooring:
      name = mooring
      symbol = "triangle"
    elif station and mooring:
      name = station, " ", mooring
      symbol = "triangle"
    else:
      name = station
      symbol = "circle"
    print "        <name>{}</name>" .format(name)
    print "        <desc>", mooring, mooring_site, operation, depth, "</desc>"
    print "        <sym>{}</sym>" .format(symbol)
    print "        <type>WPT</type>"
    print "    </wpt>"

print "</gpx>"