/******************************************************************************
 * src/wxClickText.cpp
 *
 * Implements a wxStaticText label which allows a tooltip and click events.
 *
 * based on http://forums.wxwidgets.org/viewtopic.php?t=2416&p=9955
 *
 ******************************************************************************
 * Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#include "wxClickText.h"

wxClickText::wxClickText(wxWindow *parent, wxWindowID id, const wxString &label,
                         const wxPoint& pos,
                         const wxSize& size,
                         int style, const wxString& name)
    : wxStaticText(parent, id, label, pos, size, style | wxPOPUP_WINDOW, name)
{
}

wxClickText::~wxClickText()
{
}

void wxClickText::OnMouseLeftDownEvent(wxMouseEvent& event)
{
    wxCommandEvent myevent(wxEVT_COMMAND_BUTTON_CLICKED, GetId());
    wxPostEvent(this, myevent);
    event.Skip();
}

BEGIN_EVENT_TABLE(wxClickText, wxStaticText)
    EVT_LEFT_DOWN(wxClickText::OnMouseLeftDownEvent)
END_EVENT_TABLE()
