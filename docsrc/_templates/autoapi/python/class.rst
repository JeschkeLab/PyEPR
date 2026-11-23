{# Class template: groups members under "Attributes" and "Methods" headings
   instead of listing them in source order. Properties are treated as
   attributes. Members stay nested in the py:class directive so their
   signatures and cross-reference targets are unchanged; reST does not allow
   section titles inside a directive, so the headings are rubrics.
   Only members defined (or overridden) by this class are documented in full;
   inherited members are listed as links to the class that defines them. #}
{% macro inherited_links(members, role) %}
   {% for group in members|groupby("obj.inherited_from.full_name") %}

   .. container:: autoapi-inherited

      **Inherited from** :py:class:`~{{ group.grouper }}`: {% for member in group.list %}:py:{{ role }}:`~{{ group.grouper }}.{{ member.short_name }}`{% if not loop.last %}, {% endif %}{% endfor %}

   {% endfor %}
{% endmacro %}
{% if obj.display %}
   {% if is_own_page %}
{{ obj.id }}
{{ "=" * obj.id | length }}

   {% endif %}
   {% set visible_children = obj.children|selectattr("display")|list %}
   {% set own_page_children = visible_children|selectattr("type", "in", own_page_types)|list %}
   {% set inline_children = visible_children|rejectattr("type", "in", own_page_types)|list %}
   {% set inline_attributes = inline_children|selectattr("type", "in", ["attribute", "property"])|list %}
   {% set inline_methods = inline_children|selectattr("type", "equalto", "method")|list %}
   {% set inline_other = inline_children|rejectattr("type", "in", ["attribute", "property", "method"])|list %}
   {% if is_own_page and own_page_children %}
.. toctree::
   :hidden:

      {% for child in own_page_children %}
   {{ child.include_path }}
      {% endfor %}

   {% endif %}
   {% if is_own_page %}
{# Declare the module so sphinx.ext.viewcode can add [source] links;
   without it own-page classes have no module and get no links. #}
.. py:currentmodule:: {{ obj.id[:-(obj.qual_name|length + 1)] }}

   {% endif %}
.. py:{{ obj.type }}:: {% if is_own_page %}{{ obj.qual_name }}{% else %}{{ obj.short_name }}{% endif %}{% if obj.type_params %}[{{ obj.type_params }}]{% endif %}{% if obj.args %}({{ obj.args }}){% endif %}

   {% for (args, return_annotation) in obj.overloads %}
      {{ " " * (obj.type | length) }}   {{ obj.short_name }}{% if args %}({{ args }}){% endif %}

   {% endfor %}
   {% if obj.bases %}
      {% if "show-inheritance" in autoapi_options %}

   Bases: {% for base in obj.bases %}{{ base|link_objs }}{% if not loop.last %}, {% endif %}{% endfor %}
      {% endif %}


      {% if "show-inheritance-diagram" in autoapi_options and obj.bases != ["object"] %}
   .. autoapi-inheritance-diagram:: {{ obj.obj["full_name"] }}
      :parts: 1
         {% if "private-members" in autoapi_options %}
      :private-bases:
         {% endif %}

      {% endif %}
   {% endif %}
   {% if obj.docstring %}

   {{ obj.docstring|indent(3) }}
   {% endif %}
   {% if inline_attributes %}

   .. rubric:: Attributes
      :class: autoapi-member-heading

      {% for obj_item in inline_attributes|rejectattr("inherited") %}

   {{ obj_item.render()|indent(3) }}
      {% endfor %}
{{ inherited_links(inline_attributes|selectattr("inherited")|list, "attr") }}
   {% endif %}
   {% if inline_methods %}

   .. rubric:: Methods
      :class: autoapi-member-heading

      {% for obj_item in inline_methods|rejectattr("inherited") %}

   {{ obj_item.render()|indent(3) }}
      {% endfor %}
{{ inherited_links(inline_methods|selectattr("inherited")|list, "meth") }}
   {% endif %}
   {% for obj_item in inline_other %}

   {{ obj_item.render()|indent(3) }}
   {% endfor %}
   {% if is_own_page and own_page_children %}
      {% set visible_attributes = own_page_children|selectattr("type", "in", ["attribute", "property"])|list %}
      {% if visible_attributes %}
Attributes
----------

.. autoapisummary::

         {% for attribute in visible_attributes %}
   {{ attribute.id }}
         {% endfor %}


      {% endif %}
      {% set visible_exceptions = own_page_children|selectattr("type", "equalto", "exception")|list %}
      {% if visible_exceptions %}
Exceptions
----------

.. autoapisummary::

         {% for exception in visible_exceptions %}
   {{ exception.id }}
         {% endfor %}


      {% endif %}
      {% set visible_classes = own_page_children|selectattr("type", "equalto", "class")|list %}
      {% if visible_classes %}
Classes
-------

.. autoapisummary::

         {% for klass in visible_classes %}
   {{ klass.id }}
         {% endfor %}


      {% endif %}
      {% set visible_methods = own_page_children|selectattr("type", "equalto", "method")|list %}
      {% if visible_methods %}
Methods
-------

.. autoapisummary::

            {% for method in visible_methods %}
   {{ method.id }}
            {% endfor %}


      {% endif %}
   {% endif %}
{% endif %}
