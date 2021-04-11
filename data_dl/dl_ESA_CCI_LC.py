
import cdsapi


try:
    c = cdsapi.Client()
    c.retrieve(
        'satellite-land-cover',
        {
            'variable': 'all',
            'format': 'tgz',
            'year': [
                '2014', '2015', '2016',
                '2017', '2018',
            ],
            'version': [
                'v2.0.7cds', 'v2.1.1',
            ],
        },
        'download.tar.gz')
except Exception as e:
    print("\n\nThe following exception was raised:\n\t%s\n\n" % e)
    print(("\n\nNOTE:   To run this file, the 'cdsapi' package "
           "must be correctly installed  and configured. "
           "For instructions, see https://pypi.org/project/cdsapi/\n\n"))
