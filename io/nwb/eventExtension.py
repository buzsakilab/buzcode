from pynwb.misc import AnnotationSeries
from pynwb.spec import NWBAttributeSpec


spec = NWBDatasetSpec('Event',
                    attribute=[
                        NWBAttributeSpec('timestamps', nan, 'float'),
                    ],
                    dtype=[
                        NWBDtypeSpec('foo', 'column for foo', 'int'),
                        NWBDtypeSpec('bar', 'a column for bar', 'float')
                    ])