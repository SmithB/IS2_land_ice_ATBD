function year=date2year(days)

year=(days-datenum('jan 1 2000'))/365.25 +2000;
