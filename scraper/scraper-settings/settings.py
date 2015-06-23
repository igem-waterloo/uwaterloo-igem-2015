# -*- coding: utf-8 -*-

# Scrapy settings for igempnp project
#
# For simplicity, this file contains only the most important settings by
# default. All the other settings are documented here:
#
#     http://doc.scrapy.org/en/latest/topics/settings.html
#

BOT_NAME = 'igempnp'

SPIDER_MODULES = ['igempnp.spiders']
NEWSPIDER_MODULE = 'igempnp.spiders'


DUPEFILTER_CLASS = 'scrapy.dupefilter.BaseDupeFilter'

ITEM_PIPELINES = {'igempnp.pipelines.CSVPipeline': 300 }

#DUPEFILTER_DEBUG = True

# Crawl responsibly by identifying yourself (and your website) on the user-agent
#USER_AGENT = 'igempnp (+http://www.yourdomain.com)'
