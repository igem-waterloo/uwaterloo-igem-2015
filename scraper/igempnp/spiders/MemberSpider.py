import scrapy
from scrapy.contrib.spiders import CrawlSpider, Rule
from scrapy.contrib.linkextractors import LinkExtractor

from igempnp.items import MemberItem, TrackItem

class MemberSpider(CrawlSpider):

    dont_filter = True

    name = "track"
    allowed_domains = ["igem.org"]
    start_urls = ["http://igem.org/Team_List?year=2007",
                  "http://igem.org/Team_List?year=2008",
                  "http://igem.org/Team_List?year=2009",
                  "http://igem.org/Team_List?year=2010",
                  "http://igem.org/Team_List?year=2011",
                  "http://igem.org/Team_List?year=2012",
                  "http://igem.org/Team_List?year=2013",
                  "http://igem.org/Team_List?year=2014"]

    rules = (Rule(LinkExtractor(allow=('Team\.cgi', )), callback='parse_team'),)

    def parse_team(self, response):
        team = response.xpath('//td/text()').extract()[1]
        region = response.xpath('//td/text()').extract()[11]
        year = response.xpath('//title/text()').extract()[0].split()[1]
        track = response.xpath('//table[@id="table_tracks"]/tr/td/text()')[0].extract()[16:-1]

        if track == "t been assigned to a track.":
            track = "Undefined"

        teamtrack = TrackItem()
        teamtrack['team_year'] = team+'_'+year
        teamtrack['track'] = track
        teamtrack['region'] = region

        yield teamtrack
