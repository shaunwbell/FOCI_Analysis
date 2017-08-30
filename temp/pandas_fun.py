# -*- coding: utf-8 -*-
"""
Created on Wed May 18 13:42:02 2016

@author: bell
"""
plt.style.use('ggplot')

fig = plt.figure()
plt.subplot(9, 1, 1)
plt.plot_date(pd.to_datetime(d_2003.index), d_2003['Daily Temp'],'r-')
plt.ylim([245, 285])
plt.xlim([datetime.datetime.toordinal(datetime.datetime.strptime('1-2-2003','%d-%m-%Y')), 
          datetime.datetime.toordinal(datetime.datetime.strptime('1-6-2003','%d-%m-%Y'))])
ax = plt.gca()
ax.legend(lines[:1], ['2003'], loc='lower right')
ax.xaxis.set_major_locator(MonthLocator())
ax.xaxis.set_major_formatter(ticker.NullFormatter())
ax.xaxis.set_minor_locator(MonthLocator(bymonth=range(1,13), bymonthday=15))
ax.tick_params(axis='both', which='minor', labelsize=12)

plt.subplot(9, 1, 2)
plt.plot_date(pd.to_datetime(d_2007.index), d_2007['Daily Temp'],'r-')
plt.ylim([245, 285])
plt.xlim([datetime.datetime.toordinal(datetime.datetime.strptime('1-2-2007','%d-%m-%Y')), 
          datetime.datetime.toordinal(datetime.datetime.strptime('1-6-2007','%d-%m-%Y'))])
ax = plt.gca()
ax.legend(lines[:1], ['2007'], loc='lower right')
ax.xaxis.set_major_locator(MonthLocator())
ax.xaxis.set_major_formatter(ticker.NullFormatter())
ax.xaxis.set_minor_locator(MonthLocator(bymonth=range(1,13), bymonthday=15))
ax.tick_params(axis='both', which='minor', labelsize=12)

plt.subplot(9, 1, 3)
plt.plot_date(pd.to_datetime(d_2008.index), d_2008['Daily Temp'],'r-')
plt.ylim([245, 285])
plt.xlim([datetime.datetime.toordinal(datetime.datetime.strptime('1-2-2008','%d-%m-%Y')), 
          datetime.datetime.toordinal(datetime.datetime.strptime('1-6-2008','%d-%m-%Y'))])
ax = plt.gca()
ax.legend(lines[:1], ['2008'], loc='lower right')
ax.xaxis.set_major_locator(MonthLocator())
ax.xaxis.set_major_formatter(ticker.NullFormatter())
ax.xaxis.set_minor_locator(MonthLocator(bymonth=range(1,13), bymonthday=15))
ax.tick_params(axis='both', which='minor', labelsize=12)

plt.subplot(9, 1, 4)
plt.plot_date(pd.to_datetime(d_2009.index), d_2009['Daily Temp'],'r-')
plt.ylim([245, 285])
plt.xlim([datetime.datetime.toordinal(datetime.datetime.strptime('1-2-2009','%d-%m-%Y')), 
          datetime.datetime.toordinal(datetime.datetime.strptime('1-6-2009','%d-%m-%Y'))])
ax = plt.gca()
ax.legend(lines[:1], ['2009'], loc='lower right')
ax.xaxis.set_major_locator(MonthLocator())
ax.xaxis.set_major_formatter(ticker.NullFormatter())
ax.xaxis.set_minor_locator(MonthLocator(bymonth=range(1,13), bymonthday=15))
ax.tick_params(axis='both', which='minor', labelsize=12)

plt.subplot(9, 1, 5)
plt.plot_date(pd.to_datetime(d_2010.index), d_2010['Daily Temp'],'r-')
plt.ylim([245, 285])
plt.xlim([datetime.datetime.toordinal(datetime.datetime.strptime('1-2-2010','%d-%m-%Y')), 
          datetime.datetime.toordinal(datetime.datetime.strptime('1-6-2010','%d-%m-%Y'))])
ax = plt.gca()
ax.legend(lines[:1], ['2010'], loc='lower right')
ax.xaxis.set_major_locator(MonthLocator())
ax.xaxis.set_major_formatter(ticker.NullFormatter())
ax.xaxis.set_minor_locator(MonthLocator(bymonth=range(1,13), bymonthday=15))
ax.tick_params(axis='both', which='minor', labelsize=12)

plt.subplot(9, 1, 6)
plt.plot_date(pd.to_datetime(d_2012.index), d_2012['Daily Temp'],'r-')
plt.ylim([245, 285])
plt.xlim([datetime.datetime.toordinal(datetime.datetime.strptime('1-2-2012','%d-%m-%Y')), 
          datetime.datetime.toordinal(datetime.datetime.strptime('1-6-2012','%d-%m-%Y'))])
ax = plt.gca()
ax.legend(lines[:1], ['202'], loc='lower right')
ax.xaxis.set_major_locator(MonthLocator())
ax.xaxis.set_major_formatter(ticker.NullFormatter())
ax.xaxis.set_minor_locator(MonthLocator(bymonth=range(1,13), bymonthday=15))
ax.tick_params(axis='both', which='minor', labelsize=12)

plt.subplot(9, 1, 7)
plt.plot_date(pd.to_datetime(d_2013.index), d_2013['Daily Temp'],'r-')
plt.ylim([245, 285])
plt.xlim([datetime.datetime.toordinal(datetime.datetime.strptime('1-2-2013','%d-%m-%Y')), 
          datetime.datetime.toordinal(datetime.datetime.strptime('1-6-2013','%d-%m-%Y'))])
ax = plt.gca()
ax.legend(lines[:1], ['2013'], loc='lower right')
ax.xaxis.set_major_locator(MonthLocator())
ax.xaxis.set_major_formatter(ticker.NullFormatter())
ax.xaxis.set_minor_locator(MonthLocator(bymonth=range(1,13), bymonthday=15))
ax.tick_params(axis='both', which='minor', labelsize=12)

plt.subplot(9, 1, 8)
plt.plot_date(pd.to_datetime(d_2014.index), d_2014['Daily Temp'],'r-')
plt.ylim([245, 285])
plt.xlim([datetime.datetime.toordinal(datetime.datetime.strptime('1-2-2014','%d-%m-%Y')), 
          datetime.datetime.toordinal(datetime.datetime.strptime('1-6-2014','%d-%m-%Y'))])
ax = plt.gca()
ax.legend(lines[:1], ['2014'], loc='lower right')
ax.xaxis.set_major_locator(MonthLocator())
ax.xaxis.set_major_formatter(ticker.NullFormatter())
ax.xaxis.set_minor_locator(MonthLocator(bymonth=range(1,13), bymonthday=15))
ax.tick_params(axis='both', which='minor', labelsize=12)

plt.subplot(9, 1, 9)
plt.plot_date(pd.to_datetime(d_2015.index), d_2015['Daily Temp'],'r-')
plt.ylim([245, 285])
plt.xlim([datetime.datetime.toordinal(datetime.datetime.strptime('1-2-2015','%d-%m-%Y')), 
          datetime.datetime.toordinal(datetime.datetime.strptime('1-6-2015','%d-%m-%Y'))])
ax = plt.gca()
ax.legend(lines[:1], ['2015'], loc='lower right')
ax.xaxis.set_major_locator(MonthLocator())
ax.xaxis.set_major_formatter(ticker.NullFormatter())
ax.xaxis.set_minor_locator(MonthLocator(bymonth=range(1,13), bymonthday=15))
ax.xaxis.set_minor_formatter(DateFormatter('%b'))
ax.tick_params(axis='both', which='minor', labelsize=12)