import cdsapi
import datetime as dt
c = cdsapi.Client()
def cams_solar_radiation_ts(latitude, longitude, start = None, end = None, altitude = None, filename = None):
  
  if start is None:
    start = "2005-01-01"
    print('The argument "start" is missing, default is '+ start)
  else: 
    start = dt.datetime.strptime(start, "%Y-%m-%d").strftime("%Y-%m-%d")

  if end is None:
    end = dt.date.today().strftime("%Y-%m-%d")
    print('The argument "end" is missing, default is '+ end)
  else: 
    end = dt.datetime.strptime(end, "%Y-%m-%d").strftime("%Y-%m-%d")
     
  if altitude is None:
    altitude = '-999.'
    print('The argument "altitude" is missing, default is '+ altitude)
  else: 
    altitude = str(altitude)
    
  if filename is None:
    filename = 'cams-'+str(latitude)+'_'+str(longitude)+'_'+start+'_'+end+'.csv'
    print('The argument "filename" is missing, default is '+ filename)
  else: 
    filename = str(filename)
  c.retrieve(
      'cams-solar-radiation-timeseries',
      {
          'sky_type': 'observed_cloud',
          'location': {
              'latitude': latitude,
              'longitude': longitude,
          },
          'altitude': altitude,
          'date': start+'/'+end,
          'time_step': '1day',
          'time_reference': 'universal_time',
          'format': 'csv',
      }, filename)

# cams_solar_radiation_ts(41, 15, start = '2005-01-01', end = '2023-09-01')


    
