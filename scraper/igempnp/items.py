# -*- coding: utf-8 -*-

# Define here the models for your scraped items
#
# See documentation in:
# http://doc.scrapy.org/en/latest/topics/items.html

import scrapy


class MemberItem(scrapy.Item):
    # define the fields for your item here like:
    # name = scrapy.Field()

    
    team_year = scrapy.Field()
    

class TrackItem(scrapy.Item):
    # define the fields for your item here like:
    # name = scrapy.Field()

    team_year = scrapy.Field()
    track = scrapy.Field()
    region = scrapy.Field()
    member = scrapy.Field()
    role = scrapy.Field()