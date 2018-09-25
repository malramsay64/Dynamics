#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.

"""Code to assist in the creation of figures."""

import altair as alt


def my_theme():
    return {"config": {"view": {"height": 400, "width": 600}, "axis": {"grid": False}}}


def use_my_theme():
    # register and enable the theme
    alt.themes.register("my_theme", my_theme)
    alt.themes.enable("my_theme")
